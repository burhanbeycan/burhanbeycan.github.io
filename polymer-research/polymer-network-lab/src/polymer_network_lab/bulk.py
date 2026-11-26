"""Run Polyply packing (optional) and GROMACS stages from supplied validated inputs."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys


def file_hash(path):
    digest=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            digest.update(chunk)
    return digest.hexdigest()


def identifier(value):
    if not isinstance(value,str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,60}",value):
        raise ValueError("Use a short replica or stage ID containing letters, digits, hyphens or underscores.")
    return value


def input_file(source, name):
    if not isinstance(name,str) or not name or Path(name).is_absolute():
        raise ValueError("Input file paths must be relative to their replica input folder.")
    path=(source/name).resolve()
    if not path.is_relative_to(source) or not path.is_file():
        raise ValueError(f"Missing input or file outside the replica folder: {name}")
    return path.relative_to(source)


def build_plan(config, out, gmx="gmx", polyply="polyply"):
    config=Path(config).resolve()
    out=Path(out).resolve()
    spec=json.loads(config.read_text(encoding="utf-8"))
    if not isinstance(spec,dict) or not isinstance(spec.get("replicas"),list) or not spec["replicas"]:
        raise ValueError("Define at least one replica.")
    plan={"status":"planned","config":str(config),"config_sha256":file_hash(config),"replicas":[]}
    used=set()
    for r in spec["replicas"]:
        if not isinstance(r,dict) or not isinstance(r.get("system_dir"),str) or not r["system_dir"]:
            raise ValueError("Supply a system_dir for each replica.")
        rid=identifier(r.get("id"))
        if rid in used:
            raise ValueError("Replica IDs must be unique.")
        used.add(rid)
        source=(config.parent/r["system_dir"]).resolve()
        if not source.is_dir() or out.is_relative_to(source):
            raise ValueError("Use a dedicated input folder and place the run output outside every input folder.")
        meta=r.get("evidence",spec.get("evidence"))
        if not isinstance(meta,str) or not meta.strip():
            raise ValueError("Describe the force field, validation and preparation/seed basis in evidence.")
        stages=r.get("stages",spec.get("stages"))
        if not isinstance(stages,list) or not stages:
            raise ValueError("Specify one or more stage definitions and MDP files.")
        work=out/rid
        snapshot=work/"inputs"
        def inp(name):
            return snapshot/input_file(source,name)
        top=inp(r.get("topology"))
        moves=[]
        def add(label,args,expected):
            moves.append({"step":label,"argv":[str(x) for x in args],"cwd":str(snapshot),
                          "expected":[str(x) for x in expected],"status":"planned"})
        if ("coordinates" in r)==("polyply_box_nm" in r):
            raise ValueError("Choose either prepacked coordinates or a Polyply box for each replica.")
        if "coordinates" in r:
            coords=inp(r["coordinates"])
        else:
            box=r["polyply_box_nm"]
            if not isinstance(box,list) or len(box)!=3 or any(type(x) not in (int,float) or not math.isfinite(x) or x<=0 for x in box):
                raise ValueError("Polyply box needs three positive lengths in nm.")
            coords=work/"packed.gro"
            add("pack",[polyply,"gen_coords","-p",top,"-o",coords,"-name",rid,"-box",*box],[coords])
        previous=None
        stage_names=set()
        for stage in stages:
            if not isinstance(stage,dict):
                raise ValueError("Each stage must have an ID and MDP.")
            sid=identifier(stage.get("id"))
            if sid in stage_names:
                raise ValueError("Stage IDs must be unique within the replica.")
            stage_names.add(sid)
            prefix=work/("stage-"+sid)
            tpr=prefix.with_suffix(".tpr")
            processed=work/(sid+"-processed.top")
            resolved=work/(sid+"-processed.mdp")
            args=[gmx,"grompp","-f",inp(stage.get("mdp")),"-p",top,"-c",coords,
                  "-o",tpr,"-pp",processed,"-po",resolved]
            if "index" in r:
                args.extend(["-n",inp(r["index"])])
            if "restraint_reference" in stage:
                args.extend(["-r",inp(stage["restraint_reference"])])
            use_cpt=stage.get("continue_from_previous_checkpoint",False)
            if type(use_cpt) is not bool or (use_cpt and previous is None):
                raise ValueError("Checkpoint continuation requires an earlier stage and a boolean flag.")
            if use_cpt:
                args.extend(["-t",previous.with_suffix(".cpt")])
            add("prepare-"+sid,args,[tpr,processed,resolved])
            coords=prefix.with_suffix(".gro")
            add("run-"+sid,[gmx,"mdrun","-s",tpr,"-deffnm",prefix],[coords])
            previous=prefix
        files={}
        for path in sorted(source.rglob("*")):
            if path.is_symlink():
                raise ValueError("Bundle local force-field includes as regular files inside system_dir.")
            if path.is_file():
                files[str(path.relative_to(source))]=file_hash(path)
        plan["replicas"].append({"id":rid,"source":str(source),"inputs_sha256":files,"evidence":meta,
                                 "steps":moves,"status":"planned"})
    return plan


def save_plan(out,plan):
    temp=out/"manifest.tmp"
    temp.write_text(json.dumps(plan,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    temp.replace(out/"manifest.json")


def execute(plan,out):
    out=Path(out).resolve()
    # A new directory is acquired atomically, so concurrent/repeated submissions cannot clobber runs.
    out.mkdir(parents=True,exist_ok=False)
    plan["status"]="running"
    save_plan(out,plan)
    for r in plan["replicas"]:
        loc=out/r["id"]
        step=None
        try:
            shutil.copytree(r["source"],loc/"inputs")
            actual={str(p.relative_to(loc/"inputs")):file_hash(p) for p in (loc/"inputs").rglob("*") if p.is_file()}
            if actual!=r["inputs_sha256"]:
                raise ValueError("Input files changed after planning; start a new run with the intended inputs.")
            r["status"]="running"
            for step in r["steps"]:
                exe=shutil.which(step["argv"][0])
                if not exe:
                    raise FileNotFoundError(f"Executable unavailable: {step['argv'][0]}")
                step["executable"]=str(Path(exe).resolve())
                step["executable_sha256"]=file_hash(Path(exe))
                step["status"]="running"
                save_plan(out,plan)
                with (loc/(step["step"]+".log")).open("w",encoding="utf-8") as log:
                    completed=subprocess.run([exe,*step["argv"][1:]],cwd=step["cwd"],stdout=log,stderr=subprocess.STDOUT)
                step["returncode"]=completed.returncode
                if completed.returncode:
                    raise RuntimeError(f"{step['step']} exited {completed.returncode}; see its log.")
                for target in map(Path,step["expected"]):
                    if not target.is_file() or target.stat().st_size==0:
                        raise RuntimeError(f"Expected nonempty output missing: {target}")
                step["output_sha256"]={str(Path(p).name):file_hash(Path(p)) for p in step["expected"]}
                step["status"]="completed"
            r["status"]="completed"
        except (OSError,RuntimeError,ValueError,KeyboardInterrupt) as exc:
            plan["status"]=r["status"]="failed"
            if step:
                step["status"]="failed"
            r["error"]=str(exc) or type(exc).__name__
            save_plan(out,plan)
            raise
        save_plan(out,plan)
    plan["status"]="completed"
    save_plan(out,plan)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config",type=Path)
    parser.add_argument("output",type=Path,help="New output directory.")
    parser.add_argument("--gmx",default="gmx",help="Executable name in PATH or absolute path.")
    parser.add_argument("--polyply",default="polyply",help="Executable name in PATH or absolute path.")
    parser.add_argument("--dry-run",action="store_true",help="Print plan; no tools or output directory required.")
    args=parser.parse_args()
    try:
        plan=build_plan(args.config,args.output,args.gmx,args.polyply)
        if args.dry_run:
            print(json.dumps(plan,indent=2,ensure_ascii=False))
        else:
            execute(plan,args.output)
            print(args.output.resolve()/"manifest.json")
    except (OSError,ValueError,RuntimeError) as exc:
        print(str(exc),file=sys.stderr)
        return 1
    return 0


if __name__=="__main__":
    sys.exit(main())
