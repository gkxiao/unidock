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
<p>Figure 1. Left:before the revision;  Middle: Meeko;  Right: after the revision</p>

<h2>ligand preparation</h2>
<h3>Translate SDF into Unidock-style SDF</h3>
<pre lang="shell">
    sdf2unidocksdf.py actives_final.sdf actives_prep.sdf
</pre>
<h3>Split sdf file</h3>
<pre lang="shell">
    obabel -isdf actives_prep.sdf -osdf -O actives/actives-.sdf -m
</pre>
<h3>Generate ligand index</h3>
<pre lang="shell">
    ls actives/*.sdf >> actives.index
</pre>

<h2>Docking with fast mode</h2>
<pre lang="shell">
    unidock --config dock.conf --ligand_index actives.index --dir actives_out --search_mode fast
</pre>

<h2>Post-docking</h2>
<p>Collect the docking results:</p>
<pre lang="shell">
    for file in actives_out/*.sdf
    do
      cat $file >> actives_dock.sdf
    done
</pre>
<p>Because the bond length and oriention of non-polar hydrogens are meaningless, they need to be removed and then added.</p>
<p>Remove hydrogen with openbabel:</p>
<pre lang="shell">
    obabel -isdf actives_dock.sdf -osdf -O actives_dock_noH.sdf -d
</pre>
<p>Add hydrogens with openbabel:</p>
<pre lang="shell">
    obabel -isdf actives_dock_noH.sdf -osdf -O actives_dock_addH.sdf -h
</pre>
