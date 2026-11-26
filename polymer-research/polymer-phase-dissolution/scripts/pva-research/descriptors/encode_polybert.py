"""Generate academic comparison fingerprints using local, verified polyBERT model files."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
MODEL_REPO = 'xushijie/polyBERT'
REVISION = 'e7dce434fb3eff37905dc114008660e5479ca9a8'
REQUIRED = {
    'model.safetensors': '21cfab4a3b0599852fb8ef529a67e41632246b9256ed834938ba9b46cb24148f',
    'config.json': 'ad391d9de8f9eff355554a7bb4635714250e1079f51b6b973f46d8b008be559a',
    'tokenizer.json': '89891f8a97ba93608626d704722a8243fb64c700f4e9b09edcebe794680d3c78',
    '1_Pooling/config.json': '683f6c40a16773ff51acae0782c2feadf657d547d83839cb3ceae12559674ac5'
}


def check_files(folder, download):
    from urllib.request import urlopen
    for name, expected in REQUIRED.items():
        file = folder / name
        if not file.exists() and download:
            with urlopen(f'https://huggingface.co/{MODEL_REPO}/resolve/{REVISION}/{name}', timeout=120) as response:
                content = response.read()
            if hashlib.sha256(content).hexdigest() != expected:
                raise ValueError(f'Unexpected download SHA-256 for {name}')
            file.parent.mkdir(parents=True, exist_ok=True)
            file.write_bytes(content)
        if not file.is_file():
            raise FileNotFoundError(f'Missing {file}. Supply the pinned files or --download.')
        if hashlib.sha256(file.read_bytes()).hexdigest() != expected:
            raise ValueError(f'Unexpected file SHA-256 for {name}')


def encode(folder, smiles):
    import torch
    from transformers import AutoModel
    from tokenizers import Tokenizer
    tokenizer = Tokenizer.from_file(str(folder/'tokenizer.json'))
    model = AutoModel.from_pretrained(str(folder), local_files_only=True, use_safetensors=True).eval()
    torch.set_num_threads(2)
    ids = [tokenizer.encode(s).ids for s in smiles]
    if any(len(a) > 512 for a in ids):
        raise ValueError('One repeat exceeds the model context.')
    if any(tokenizer.token_to_id('[UNK]') in a for a in ids):
        raise ValueError('One input contains an unknown model token.')
    pad = tokenizer.token_to_id('[PAD]')
    input_ids = torch.nn.utils.rnn.pad_sequence([torch.tensor(a) for a in ids], batch_first=True, padding_value=pad)
    mask = (input_ids != pad).long()
    with torch.inference_mode():
        hidden = model(input_ids=input_ids, attention_mask=mask).last_hidden_state
        pooled = (hidden*mask.unsqueeze(-1)).sum(dim=1)/mask.sum(dim=1, keepdim=True)
    values = pooled.cpu().numpy()
    if values.shape != (len(smiles), 600) or not np.isfinite(values).all():
        raise ValueError('Unexpected fingerprints.')
    return values


def normalized_input(psmiles, catalog):
    # Reuse the registered spelling for graph-identical motifs to preserve
    # comparability with frozen benchmark embeddings.
    from rdkit import Chem
    from structure import repeat_descriptors
    def key(text):
        repeat_descriptors(text)
        mol = Chem.MolFromSmiles(text)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    graph = key(psmiles)
    for row in catalog.values():
        if key(row['psmiles']) == graph:
            return row['psmiles']
    return graph.replace('*', '[*]')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model-dir', type=Path, required=True)
    p.add_argument('--download', action='store_true')
    p.add_argument('--out', type=Path, help='Output .npz path; required for custom inputs')
    p.add_argument('--repeat', help='One linear PSMILES repeat; output key is QUERY')
    args = p.parse_args()
    check_files(args.model_dir, args.download)
    catalog = json.loads((HERE/'repeat_catalog.json').read_text())
    if args.repeat:
        if args.out is None:
            p.error('--repeat requires a distinct --out path')
        if args.out.resolve() == (HERE/'polybert_embeddings.npz').resolve():
            p.error('Choose a separate output for the query')
        inputs = {'QUERY': normalized_input(args.repeat, catalog)}
    else:
        inputs = {key: row['psmiles'] for key,row in catalog.items()}
    values = encode(args.model_dir, list(inputs.values()))
    np.savez_compressed(args.out or HERE/'polybert_embeddings.npz', **dict(zip(inputs, values)))
    print(f'{len(inputs)} repeat motifs; {values.shape[1]} frozen features per motif.')


if __name__ == '__main__':
    main()
