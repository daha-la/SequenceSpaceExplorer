"""ESM-C sequence embedder. Works on any datafile (uses the Sequence column)."""

from .base import Embedder, EmbedContext, pool, to_numpy_fp32, resolve_device
from ..common import COL_ID, COL_SEQ


class EsmcEmbedder(Embedder):
    name = "esmc"
    requires_structure = False

    def tag(self, args) -> str:
        size = args.esmc_model.replace("esmc_", "")        # esmc_600m -> 600m
        return f"esmc{size}_{args.pooling}"

    def prepare(self, df, ctx: EmbedContext, args):
        # Every datafile row has a usable sequence (creation guaranteed it), so
        # ESM-C never skips; a failure later is a genuine runtime error.
        entries = [{"id": r[COL_ID], "seq": r[COL_SEQ]}
                   for _, r in df.iterrows()]
        return entries, {}

    def load_model(self, args):
        # NB: ESM-C runs best on CUDA; CPU/MPS support depends on the installed
        # `esm` build. If load fails on a non-CUDA device, that is the SDK.
        from esm.models.esmc import ESMC
        from esm.sdk.api import ESMProtein, LogitsConfig
        device = resolve_device(args.device)
        print(f"  device: {device} | model: {args.esmc_model} | pooling: {args.pooling}")
        client = ESMC.from_pretrained(args.esmc_model).to(device).eval()
        cfg = LogitsConfig(sequence=True, return_embeddings=True)
        return {"client": client, "cfg": cfg, "ESMProtein": ESMProtein,
                "device": device}

    def encode_batch(self, mc, batch, args):
        client, cfg, ESMProtein = mc["client"], mc["cfg"], mc["ESMProtein"]
        proteins = [ESMProtein(sequence=e["seq"]) for e in batch]
        try:
            pt = client.encode(proteins)
            batched = True
        except Exception:
            pt = [client.encode(p) for p in proteins]
            batched = False
        try:
            out = (client.logits(pt, cfg) if batched
                   else [client.logits(x, cfg) for x in pt])
        except Exception:
            seq = pt if isinstance(pt, list) else [pt]
            out = [client.logits(x, cfg) for x in seq]
        if not isinstance(out, list):
            out = list(out)
        return [to_numpy_fp32(pool(o.embeddings[-1], args.pooling)) for o in out]
