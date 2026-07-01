"""SaProt structure embedder (Foldseek-only). Embeds paired AA + 3Di tokens."""

from .base import Embedder, EmbedContext, pool, to_numpy_fp32, resolve_device, use_fp16
from ..common import COL_ID, COL_SEQ
from . import structure_input


def _saprot_seq(aa: str, di: str):
    """Each residue = AA (upper) + 3Di (lower), e.g. 'MdAaKp'. None on mismatch."""
    if len(aa) != len(di):
        return None
    return "".join(f"{a}{d.lower()}" for a, d in zip(aa, di))


class SaProtEmbedder(Embedder):
    name = "saprot"
    requires_structure = True

    def tag(self, args) -> str:
        return f"saprot_{args.pooling}"

    def prepare(self, df, ctx: EmbedContext, args):
        di_by_id, skip = structure_input.compute_3di(df, ctx, args.max_residues)
        skip["len_mismatch"] = []
        aa_by_id = {r[COL_ID]: r[COL_SEQ] for _, r in df.iterrows()}
        entries = []
        for rid, di in di_by_id.items():
            sp = _saprot_seq(aa_by_id.get(rid, ""), di)
            if sp is None:
                skip["len_mismatch"].append(rid)
                continue
            entries.append({"id": rid, "saprot_seq": sp})
        return entries, skip

    def load_model(self, args):
        from transformers import EsmTokenizer, EsmModel
        device = resolve_device(args.device)
        print(f"  device: {device} | SaProt ({args.saprot_checkpoint})")
        tok = EsmTokenizer.from_pretrained(args.saprot_checkpoint)
        model = EsmModel.from_pretrained(args.saprot_checkpoint).to(device).eval()
        if use_fp16(device):
            model = model.half()
        return {"tok": tok, "model": model, "device": device}

    def encode_batch(self, mc, batch, args):
        tok, model, device = mc["tok"], mc["model"], mc["device"]
        seqs = [e["saprot_seq"] for e in batch]
        enc = tok(seqs, padding=True, return_tensors="pt").to(device)
        out = model(**enc)
        hidden = out.last_hidden_state
        vecs = []
        for i in range(len(batch)):
            valid = enc["attention_mask"][i].bool().clone()
            idxs = valid.nonzero(as_tuple=True)[0]
            if len(idxs) >= 2:
                valid[idxs[0]] = False       # <cls>
                valid[idxs[-1]] = False       # <eos>
            vecs.append(to_numpy_fp32(pool(hidden[i], args.pooling, valid)))
        return vecs
