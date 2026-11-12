// Compiled by `tsc -b` (npm run build) and never imported: a fact without its
// source does not type-check. Each line below must fail to compile; if the type
// ever lets one through, its @ts-expect-error is unused and the build stops.

import { chosen, measured, notRecorded, recorded, type Fact } from './facts.ts'

// @ts-expect-error a recorded fact without its place
export const noPlace: Fact = { kind: 'recorded', value: 'scanvi' }
// @ts-expect-error a recorded null without where the converter looked
export const noSearch: Fact = { kind: 'recorded', value: null }
// @ts-expect-error a measured fact without the matrix it was measured on
export const noOf: Fact = { kind: 'measured', value: true, how: 'sampled' }
// @ts-expect-error a measured fact without how
export const noHow: Fact = { kind: 'measured', value: true, of: 'X' }
// @ts-expect-error a chosen fact without who chose
export const noBy: Fact = { kind: 'chosen', value: 'cell_type' }
// @ts-expect-error a kind the table does not have
export const inferred: Fact = { kind: 'inferred', value: 'harmony', at: 'obsm/X_harmony' }
// @ts-expect-error the constructor needs the place
export const noPlace2 = recorded('scanvi')
// @ts-expect-error and where the converter looked
export const noSearch2 = notRecorded()
// @ts-expect-error and of what, and how
export const noHow2 = measured(true, 'X')
// @ts-expect-error and by whom
export const noBy2 = chosen('cell_type')
