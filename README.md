<h2>unidock_ligand_preparation</h2>
<p>build_torsion_tree.py is a tool in <a href="https://github.com/dptech-corp/Uni-Dock/tree/mcdock/unidock_tools">unidock_tools</a> used to prepare the unidock sdf.Here we compare the features before and after the modification of <a href='https://github.com/dptech-corp/Uni-Dock/pull/58'>#58</a>.</p>

<blockquote cite="https://www.huxley.net/bnw/four.html](https://github.com/dptech-corp/Uni-Dock/pull/58">
<p>Update ligand torsion tree construction rules</p>

<p>This commit improves the definition of rotatable bonds in the ligand torsion tree construction. Changes include:</p>
<ul>
    <li>Excluding terminal methyl, hydroxyl, and amine groups from being defined as rotatable single bonds.</li>
    <li>Preventing single bonds within peptide bond planes and adjacent conjugated systems from being classified as rotatable.</li>
</ul>
<p>These updates address issues where inappropriate torsional degrees of freedom were assigned, leading to unrealistic ligand conformations during the docking process.</p>
</blockquote>
<p>Here is an example to show the difference:</p>

![tree comparison](https://github.com/gkxiao/unidock_ligand_preparation/blob/main/build_tree_58_rev.jpg)
<p>Figure. Left:before the revision; Middle: Meeko; Right: after the revision</p>
