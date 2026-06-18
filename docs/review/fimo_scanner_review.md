# FIMO Scanner Review

The scanner should be described as a **FIMO motif-match scanner** and the plot as a **FIMO motif-match profile**. It should not be called a structural binding-energy profile, ModCRE affinity profile, or TF occupancy profile.

## Current Behavior

- Reads selected usable PWM matrices from SQLite.
- Writes a temporary combined MEME file.
- Writes the pasted sequence as temporary FASTA.
- Calls FIMO with an argument list, not a shell string.
- Parses `fimo.tsv` and reports motif, source, start, end, strand, matched sequence, score, p-value, and q-value.
- Builds a per-position profile using best `-log10(p-value)` and overlapping hit counts.
- Uses a uniform DNA background: `A 0.25 C 0.25 G 0.25 T 0.25`.

## Review Notes

- Uniform background is acceptable for V1 but must be documented because GC/AT-rich inputs can shift FIMO p-values.
- Public deployment needs request-size limits, rate limits, timeouts, and ideally worker isolation/queueing.
- Downloaded CSVs are summaries of FIMO motif matches, not experimental binding proof.
