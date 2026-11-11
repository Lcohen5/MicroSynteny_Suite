# species_taxonomy.tsv — classification sources

Background for `visualization/species_taxonomy.tsv`. Read this before changing a
classification; several of the assignments are contested in the literature and
the choices below are deliberate.

---

## Columns

| Column | Meaning |
|---|---|
| `species` | Species name as it appears in the `species` column of a v7 TSV. Must match exactly. |
| `abbrev` | Four-letter short form, for compact axis labels. |
| `subgenus` | *Sophophora*, *Drosophila*, or *Dorsilopha*. |
| `group` | Species group. |
| `subgroup` | Species subgroup. Where no finer division is applied, repeats the group name. |
| `order` | Display order, top to bottom. Replaces the former hardcoded `SPECIES_ORDER` vector. |
| `notes` | Free text. Not parsed. |

A level with fewer than two species present in a given dataset produces no plot —
there is nothing to draw a ribbon between. Those skips should be logged, not silent.

---

## `order` is empirical, not taxonomic

The `order` column reproduces, position for position, the tip order of the GEP
phylogenetic tree derived from 36 RefSeq *Drosophila* genomes. It was verified
against the tree and against the `SPECIES_ORDER` vector previously hardcoded in
`comparative_synteny_plot.R`; all three agree at all 36 positions.

**This matters for interpretation.** Ribbons connect species that are *adjacent
in this order*, so `order` is not merely cosmetic — changing it changes which
orthology relationships the figure asserts. Do not reorder to match a published
taxonomy without re-rendering every figure.

Note in particular that *D. ficusphila* sits at position 9, immediately outside
the *melanogaster* subgroup and ahead of the other Oriental-lineage subgroups.
Several published topologies instead place *ficusphila* with *elegans* and
*rhopaloa*. Position 9 is what the GEP tree recovers and is retained on that
basis.

---

## Classification convention followed

**The *montium* and *ananassae* lineages are treated as species groups, not as
subgroups of the *melanogaster* group.**

This follows Da Lage et al. (2007), who proposed restricting the *melanogaster*
group to the *melanogaster* subgroup plus the Oriental subgroups and raising
*ananassae* and *montium* to group rank, and Yassin (2013, 2018), who formally
subdivided the *montium* group into seven subgroups: *kikkawai*, *montium*,
*orosa*, *parvula*, *punjabiensis*, *serrata*, and *seguyi*.

This classification "has been adopted by *Drosophila* taxonomic databases
(Bächli 2020; Toda 2020) and several recent comparative analyses"
(Conner et al. 2021).

**It is not universal.** Suvorov et al. (2019) describe the elevation as
"followed by some authors … but has not been widely adopted," and much of the
older literature students will encounter treats *montium* and *ananassae* as
subgroups within a broadly defined *melanogaster* group. If you prefer the
traditional treatment, set `group = melanogaster` for rows 17–20 and
`subgroup = montium` / `ananassae` respectively; nothing else changes.

Consequence of the convention used here: *montium* and *ananassae* become
plottable groups (2 species each), and *kikkawai* and *serrata* fall into
*different* subgroups, so no *montium* subgroup plot is produced.

---

## Per-assignment sources

**Subgenera.** *D. busckii* is placed in subgenus *Dorsilopha*, positioned
ancestral to subgenus *Drosophila*. Its `group`/`subgroup` values are
placeholders — *Dorsilopha* is represented by a single species here.

**obscura group.** Six subgroups are recognised: *affinis*, *microlabis*,
*obscura*, *pseudoobscura*, *subobscura*, *sinobscura*. *D. guanche* and
*D. subobscura* are both in the ***subobscura*** subgroup, together with
*D. madeirensis* — not the *obscura* subgroup, which is a common slip given the
shared name. *D. obscura* itself is the type species of the *obscura* subgroup.

**Oriental lineage.** The *elegans*, *ficusphila*, *eugracilis*, *takahashii*,
*rhopaloa* and *favohirta* subgroups are all within the *melanogaster* species
group. They are subgroups, not groups.

**suzukii subgroup.** *D. suzukii* and *D. biarmipes* are confirmed members,
sister to the *takahashii* subgroup. *D. subpulchrella* is assigned here as
*D. suzukii*'s sister species; this is well established but is the one
assignment in the table not backed by a source consulted directly.

**immigrans group.** Subgenus *Drosophila*, within the immigrans-tripunctata
radiation. *D. albomicans* is in the *nasuta* subgroup with *D. nasuta*.

**Hawaiian Drosophila.** A monophyletic clade nested within *Drosophila*, not a
separate subgenus. The *grimshawi* species group is a component of the
picture-wing / nudidrosophila / ateledrosophila (PNA) clade. Conventions differ
as to whether picture-wing is the group and *grimshawi* the subgroup or the
reverse; picture-wing is the more inclusive label either way, so it is placed at
`group` here.

---

## References

- Bächli G. (2020). *TaxoDros: The Database on Taxonomy of Drosophilidae.*
- Conner W.R. et al. (2021). A phylogeny for the *Drosophila montium* species
  group: a model clade for comparative analyses. *Mol Phylogenet Evol.*
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7946709/
- Da Lage J.-L. et al. (2007). A phylogeny of Drosophilidae using the *Amyrel*
  gene: questioning the *Drosophila melanogaster* species group boundaries.
  *J Zool Syst Evol Res* 45(1):47–63.
  https://onlinelibrary.wiley.com/doi/10.1111/j.1439-0469.2006.00389.x
- Suvorov A. et al. (2019). Phylogenetic position of the *Drosophila fima* and
  *dentissima* lineages, and the status of the *D. melanogaster* species group.
  *Mol Phylogenet Evol.* https://pmc.ncbi.nlm.nih.gov/articles/PMC6708496/
- Toda M.J. (2020). *DrosWLD-Species: Drosophilidae of the World.*
- Yassin A. (2013, 2018). Revisions of the *montium* species group.
- *Drosophila* montium whole-genome resource (2020). *G3* 10(5):1443–1455.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7202002/
- Hawaiian Drosophila clade overview.
  https://en.wikipedia.org/wiki/Hawaiian_Drosophila
- *Drosophila obscura* species group subgroup structure.
  https://en.wikipedia.org/wiki/Drosophila_obscura_species_group
- *Dorsilopha.* https://en.wikipedia.org/wiki/Dorsilopha
- *Drosophila immigrans* species group.
  https://en.wikipedia.org/wiki/Drosophila_immigrans_species_group
- Ometto L. et al. (2013). Linking genomics and ecology to investigate the
  complex evolution of an invasive *Drosophila* pest. *Genome Biol Evol*
  5(4):745–757. https://academic.oup.com/gbe/article/5/4/745/587166

The Wikipedia entries above were used for subgroup membership that is stable and
uncontested; the contested rank questions are sourced to the primary literature.

---

## Species set

The 36 species here are those in the GEP 36-genome tree. The GEP annotation
protocol paper (Rele et al., *F1000Research* 11:1579) works from a 28-species
subset; the extra eight are present in the tree and in `SPECIES_ORDER` and are
retained so that curated data for them renders without editing the table.
