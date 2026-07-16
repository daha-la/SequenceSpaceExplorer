"""ProstT5 structure embedder (Foldseek-only). Embeds the 3Di sequence."""

from .base import Embedder, EmbedContext, pool, to_numpy_fp32, resolve_device, use_fp16
from . import structure_input


class ProstT5Embedder(Embedder):
    name = "prostt5"
    requires_structure = True

    def tag(self, args) -> str:
        return f"prostt5_{args.pooling}"

    def prepare(self, df, ctx: EmbedContext, args):
        di_by_id, skip = structure_input.compute_3di(df, ctx, args.max_residues)
        entries = [{"id": rid, "di_seq": di} for rid, di in di_by_id.items()]
        return entries, skip

    def load_model(self, args):
        from transformers import T5Tokenizer, T5EncoderModel
        device = resolve_device(args.device)
        print(f"  device: {device} | ProstT5 ({args.prostt5_checkpoint})")
        tok = T5Tokenizer.from_pretrained(args.prostt5_checkpoint, do_lower_case=False)
        model = T5EncoderModel.from_pretrained(args.prostt5_checkpoint).to(device).eval()
        if use_fp16(device):
            model = model.half()
        return {"tok": tok, "model": model, "device": device}

    def encode_batch(self, mc, batch, args):
        tok, model, device = mc["tok"], mc["model"], mc["device"]
        seqs = [f"<fold2AA> " + " ".join(e["di_seq"].lower()) for e in batch]
        enc = tok.batch_encode_plus(seqs, add_special_tokens=True,
                                    padding="longest", return_tensors="pt").to(device)
        out = model(input_ids=enc["input_ids"], attention_mask=enc["attention_mask"])
        hidden = out.last_hidden_state
        vecs = []
        for i in range(len(batch)):
            valid = enc["attention_mask"][i].bool().clone()
            idxs = valid.nonzero(as_tuple=True)[0]
            if len(idxs) >= 2:
                valid[idxs[0]] = False      # drop <fold2AA>
                valid[idxs[-1]] = False      # drop </s>
            vecs.append(to_numpy_fp32(pool(hidden[i], args.pooling, valid)))
        return vecs
