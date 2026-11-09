// What the command's exit status means.
//
// One number per kind of failure, so that a pipeline can tell a typo in its own
// call (2) from a file it was handed that is not an object (3), from an object
// that lacks what the call named (4), from one that cannot become a bundle at
// all (5), from a disk that would not take the result (6). The README lists
// the same table; keep them together.

export const EXIT = {
  ok: 0,
  /** A bug, or a Node that cannot run the command at all. */
  internal: 1,
  /** The call itself is malformed: an unknown flag, a missing value, a stray path. */
  usage: 2,
  /** The input is missing, not HDF5, or not an .h5ad this reader understands. */
  input: 3,
  /** A role names a column or embedding the object lacks, or one that cannot play it. */
  roles: 4,
  /** The object cannot become a bundle: scaled data only, no embedding, no annotation. */
  build: 5,
  /** The output could not be written. */
  output: 6,
} as const

export class CliError extends Error {
  exit: number
  constructor(message: string, exit: number) {
    super(message)
    this.exit = exit
  }
}
