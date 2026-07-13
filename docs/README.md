# docs/

Longer-form documentation that doesn't fit comfortably in a folder README.
Two different audiences, two different documents — don't confuse one for
the other:

| File | Audience | Contents |
|---|---|---|
| [`visualizer_guide.md`](visualizer_guide.md) | Anyone using the app | Every panel and tool in the visualizer (coordinates, filters, layers, colour, export, Boltz-2 structure prediction, RMSD), explained. |
| [`SSE_Visualizer_User_Guide.docx`](SSE_Visualizer_User_Guide.docx) | Anyone using the app | The same guide as a Word document, for anyone who'd rather read/annotate it outside a Markdown viewer (e.g. in a lab notebook or shared with a non-technical collaborator). |
| [`SSE_datafile_spec.md`](SSE_datafile_spec.md) | Anyone extending SSE | The internal engineering contract for the `.sse.tsv` format and every tool's read/write guarantees — design rationale, decided/open/deferred items included. What the `spec §X` comments throughout `sse_tools/` and `scripts/` cite. Not required reading for using SSE. |

`visualizer_guide.md` and `SSE_Visualizer_User_Guide.docx` are equivalent
in content — the Markdown is the source of truth; the Word document is
generated from it. If you're editing the guide, edit the Markdown and
regenerate the Word document from it, rather than editing the two
independently.
