"""Generate the review figure from the exact graph and aligned GUI drawing data."""
from pathlib import Path
from xml.etree import ElementTree as ET
from polymer_chemistry_lab.end_groups import petox_candidates


def node(parent, tag, attributes, label=None):
    item = ET.SubElement(parent, tag, {k:str(v) for k,v in attributes.items()})
    if label is not None:
        item.text = label
    return item


def main():
    palette = {"alpha":"#933f79", "omega":"#087d70"}
    rows = petox_candidates(3)
    points = [a for row in rows for a in row["view"]["atoms"]]
    xmin, xmax = min(a["x"] for a in points), max(a["x"] for a in points)
    ymin, ymax = min(a["y"] for a in points), max(a["y"] for a in points)
    scale = min(402 / (xmax-xmin), 119 / (ymax-ymin))
    root = ET.Element("svg", {"xmlns":"http://www.w3.org/2000/svg","viewBox":"0 0 812 558",
                              "role":"img","aria-labelledby":"title desc","font-family":"Arial,sans-serif"})
    node(root,"title",{"id":"title"},"PEtOx end groups at DP 3")
    node(root,"desc",{"id":"desc"},"Three aligned covalent structures. Magenta marks the N-terminal alpha site, green the carbon-terminal omega site. Masses refer to isolated, neutral molecules.")
    node(root,"rect",{"width":812,"height":558,"fill":"white"})
    node(root,"text",{"x":22,"y":36,"font-size":23,"fill":"#213847"},"PEtOx end groups at DP 3")
    node(root,"text",{"x":22,"y":60,"font-size":14,"fill":"#455e6b"},"Pink = α nitrogen end · Green = ω carbon end · Unlabelled skeletal vertices are carbons")
    for row_idx,row in enumerate(rows):
        top = 75 + row_idx*153
        node(root,"line",{"x1":22,"x2":790,"y1":top,"y2":top,"stroke":"#d7e1e6"})
        node(root,"text",{"x":22,"y":top+38,"font-size":19,"fill":"#243742"},row["name"])
        for j, value in enumerate([row["formula"], f'Neutral exact mass   {row["neutral_Da"]:.6f} Da',
                                   f'[M+H]⁺   {row["mz_H"]:.6f} m/z']):
            node(root,"text",{"x":22,"y":top+66+j*24,"font-size":15,"fill":"#36515e"},value)
        points = row["view"]["atoms"]
        pos = {a["index"]:(round(359+(a["x"]-xmin)*scale,1),round(top+16+(ymax-a["y"])*scale,1)) for a in points}
        for a in points:
            if "end_role" in a:
                x,y = pos[a["index"]]
                node(root,"circle",{"cx":x,"cy":y,"r":9,"fill":palette[a["end_role"]],"opacity":.20})
        for b in row["view"]["bonds"]:
            x1,y1 = pos[b["begin"]]
            x2,y2 = pos[b["end"]]
            dx,dy = x2-x1,y2-y1
            length = (dx*dx+dy*dy)**.5
            for offset in ([-1.65,1.65] if b["order"]==2 else [0]):
                px,py = -dy/length*offset,dx/length*offset
                node(root,"line",{"x1":round(x1+px,1),"y1":round(y1+py,1),"x2":round(x2+px,1),
                                   "y2":round(y2+py,1),"stroke":palette.get(b.get("end_role"),"#354952"),"stroke-width":1.7})
        for a in points:
            if a["label"]:
                x,y=pos[a["index"]]
                node(root,"text",{"x":x,"y":y+4,"text-anchor":"middle","font-size":13,
                                   "fill":palette.get(a.get("end_role"),"#213740"),"stroke":"white",
                                   "stroke-width":3,"paint-order":"stroke"},a["label"])
    node(root,"text",{"x":22,"y":548,"font-size":12,"fill":"#455e6b"},"Covalent models · prospective adduct masses · literature end assignments: Morgan et al., 2022, DOI 10.1039/D2PY00437B")
    target=Path(__file__).resolve().parents[1]/"docs"/"assets"/"petox-ends.svg"
    ET.indent(root)
    ET.ElementTree(root).write(target, encoding="utf-8", xml_declaration=True)


if __name__ == "__main__":
    main()
