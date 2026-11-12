// `scrnaseq-lab`: the commands, their help, and their exit status.

import { readFileSync } from 'node:fs'
import { parseConvertArgs } from './args.ts'
import { convert } from './convert.ts'
import { CliError, EXIT } from './errors.ts'

const USAGE = 'usage: scrnaseq-lab convert IN.h5ad OUT.zip [options]   (scrnaseq-lab --help for the options)'

const HELP = `scrnaseq-lab convert IN.h5ad OUT.zip [options]

Turns an annotated AnnData object into the one .zip scRNA-seq Studio opens.
Each part is the bundle the Lab app's own builder writes, plus "sentinels" and
the command's notes in its meta.json; collection.json holds the app's fields
plus "sentinels" and a sha256 for each part.

Roles. A flag that is given is used as given. A flag that is not given is
guessed, and the guess is printed on stderr as the flag that would give it:
  --cluster COL        the cell annotation
  --sample COL         the sample column, or none
  --condition COL      the design factor Studio groups by first, or none
  --extra COL          one more column to carry (repeatable). By default every
                       other design factor; none carries no extra column
  --interaction A,B    carry the interaction of A and B as a column of its own
                       (repeatable). By default the two carried factors'
                       interaction when they cross; none for no interaction
  --embedding KEY      the obsm embedding that opens first
  --sentinel LEVEL     a level of --cluster that is a sentinel, not a cell type
                       (repeatable). By default EXCLUDED and UNRESOLVED where
                       present; none for no sentinel
  --label TEXT         the dataset's name. By default the file name

A design factor is a column that takes one value in every sample. Every one is
carried unless a flag leaves it out, and one left out is named in the notes.

A sample whose cells hold more than one level of --condition is split into one
sample per level (sample|level), so every cell keeps its own; with no sample
column, the cells of each level are one sample. Beside the embedding that opens
first, only two-column embeddings travel, each set of coordinates once. Of
the run's records in uns, only the identity travels, as "run": tool, version,
commit, stage, run_key and job, each as a string of 256 characters or fewer.
The notes name all three. Only the input file is read: no link out of it is
followed, and the notes name each one met.

Output: one .zip holding collection.json and parts/*.zip, each part a complete
scrnaseq-studio/bundle@1. The same input and flags give the same bytes. The
file's sha256 is in the summary, and collection.json lists each part's.

stdout: one JSON summary. stderr: the guesses, the notes and the progress.

Exit status:
  0  converted
  1  internal error, or a Node that cannot run the command
  2  usage: an unknown option, a missing value, a stray path, or an output
     that is the input file (by path, symlink or hard link)
  3  input: missing, not HDF5, or not an .h5ad this reader understands
  4  roles: a flag names a column or embedding the object lacks or cannot use,
     or a --sentinel level that no cell carries
  5  the object cannot become a bundle (scaled data only, no embedding, ...)
  6  output: the file could not be written

Large objects: NODE_OPTIONS=--max-old-space-size=8192 scrnaseq-lab convert ...`

function version(): string {
  try {
    return JSON.parse(readFileSync(new URL('../../package.json', import.meta.url), 'utf8')).version
  } catch {
    return 'unknown'
  }
}

export async function main(argv: string[]): Promise<number> {
  const out = (s: string) => process.stdout.write(`${s}\n`)
  const err = (s: string) => process.stderr.write(`${s}\n`)
  const [cmd, ...rest] = argv

  if (cmd === undefined) { err(USAGE); return EXIT.usage }
  if (cmd === '--help' || cmd === '-h' || cmd === 'help') { out(HELP); return EXIT.ok }
  if (cmd === '--version') { out(version()); return EXIT.ok }
  if (cmd !== 'convert') { err(`scrnaseq-lab: unknown command "${cmd}". ${USAGE}`); return EXIT.usage }

  try {
    const args = parseConvertArgs(rest)
    if (args.help) { out(HELP); return EXIT.ok }
    out(JSON.stringify(await convert(args, err), null, 1))
    return EXIT.ok
  } catch (e) {
    if (e instanceof CliError) {
      err(`scrnaseq-lab: ${e.message}`)
      return e.exit
    }
    err(`scrnaseq-lab: internal error: ${(e as Error)?.stack ?? String(e)}`)
    return EXIT.internal
  }
}
