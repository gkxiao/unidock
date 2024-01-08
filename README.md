<h2>Introduction</h2>
<p>build_torsion_tree.py is a tool in <a href="https://github.com/dptech-corp/Uni-Dock/tree/mcdock/unidock_tools">unidock_tools</a> (mdock branch) used to prepare the unidock sdf.Here we compare the features before and after the modification of <a href='https://github.com/dptech-corp/Uni-Dock/pull/58'>#58</a>.</p>

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

<h2>Preapre unidock SDF with <a href="https://github.com/dptech-corp/Uni-Dock/tree/mcdock/unidock_tools">unidock_tools</a> (mcdock branch)</h2>
<h3>Translate SDF into Unidock-style SDF</h3>
<pre lang="shell">
sdf2unidocksdf.py actives_final.sdf actives_prep.sdf
</pre>
<h3>Split sdf file</h3>
<pre lang="shell">
obabel -isdf actives_prep.sdf -osdf -O actives/actives-.sdf -m
</pre>
<p>You can use shell command to split the SDF：</p>
<pre line="1" lang="shell">
csplit actives_prep.sdf /"\$\$\$\$"/ --suppress-matched --prefix=actives- --suffix-format=%05d.sdf '{*}'
</pre>
<p>It will generate files such as actives-00001.sdf,actives-00002.sdf....</p>
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
<h2>Preapre unidock SDF with <a href="https://github.com/dptech-corp/Uni-Dock/tree/main/unidock_tools">unidock_tools </a>(main branch) </h2>
<pre lang="python">
from unidock_tools.ligandPrepare import prepare_ligands
ligs=['4fv2_ligand.sdf',]
prepare_ligands(ligs,output_dir='.')
</pre>
<h2>Bias docking</h2>
<p>The user must prepare a bias parameter file (BPF) containing all the information for the different biases to be applied. The BPF contains one line for each bias, with the following parameters: </p>
<ul>
   <li>(x, y, z) coordinates in Å,</li> 
   <li>energy reward ( V set ) in kcal/mol,</li> 
   <li>decay radius ( r ) in Å,</li>
   <li>type of bias ( don , acc , aro or map ).</li> 
</ul>
<p>The following type of biases are available:</p>
<ul>
   <li>hydrogen bond donor = don , </li> 
   <li>hydrogen bond acceptor = acc , </li> 
   <li>aromatic = aro, </li>
   <li>specific bias according to the desired map = map </li> 
</ul>
<p>Two examples of bias parameter files are shown below: </p>
<p>i) Bias parameter file with traditional interaction biases</p> 
<pre lang="python">
x y z Vset r type 
33.682 36.327 34.471 -1.50 0.80 acc 
34.557 36.028 31.933 -2.00 0.60 don 
36.905 36.534 30.560 -1.75 1.00 aro
</pre>
<p>ii) Bias parameter file with specific map biases</p>
<pre lang="python">
x y z Vset r type 
5.100 1.785 20.019 -1.75 1.10 map 
9.459 2.075 24.527 -2.00 0.80 map
</pre>
<p>Bias parameter file format requirements</p>
<ul>
   <li>All lines must have 6 columns. The columns must be space or tab separated.</li> 
   <li>Lines are ignored if the first column is not numeric (e.g., header with titles - x, y, z, V set , r, type - in the examples above). </li> 
   <li>The first three columns define the x,y,z coordinates of the bias site center, in Å.</li>
   <li>The fourth column corresponds to the energy reward ( V set ), in kcal/mol, to be applied at the bias site center. It has to be a negative number . If there is no thermodynamic information for the site, a reasonable value is -2.0 kcal/mol, which sets a relatively strong bias.</li>
   <li>The fifth column is the radius ( r ) of the bias site, in Å. It controls the extent of energy reward through space according to a Gaussian function. A reasonable value range that builds a well defined bias site is between 0.6 and 1.2 Å.</li>
   <li>The last column indicates the type of bias and, in consequence, which energy maps will be modified: (1)acc modifies NA and OA maps;(2) don modifies HD maps; (3)aro creates an ad hoc new map (AC, aromatic center map) ;(4)map modifies the energy map specified in the -m argument.IMPORTANT: map biases
cannot be combined with other types of biases ( don , acc , aro ) in the same execution of the program.</li>
</ul>
<p>Bias-docking example:</p>
<pre lang="shell">
unidock --config dock.conf --ligand_index actives.index --dir actives_out --search_mode fast --bias hinge_ph4.bpf
</pre>
<h2>Multi-Conformation Docking (MCDOCK)</h2>
<p>Tutorial: <a href="https://nb.bohrium.dp.tech/detail/91221652314"></a>https://nb.bohrium.dp.tech/detail/91221652314</a></p>

<p>The rigid- and flexible-docking is achieved by the using fragment information.</p>
<p>This is a example of fragment information for flexible-docking：</p>
<pre lang="shell">
&gt;  &lt;fragInfo&gt;
1 2 3 4 5 30 31 32 33 34
6 7 8
9 10 11 12
13 14 15 16 21 22 23 28 29
17 18 19 20
24 25 26 27
35
36
37 38 39
40 41 42 43
44 45 46 47
</pre>

<p>This is a example of fragment information for rigid-docking：</p>
<pre lang="shell">
&gt;  &lt;fragInfo&gt;
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
</pre>
<p>If you need to relax pose, make sure to use flexible fragment information and carry out a local only search：</p>
<pre lang="shell">
unidock --config dock.conf --ligand_index pose.index --dir relax_out --exhaustiveness 512 --max_step 40 --num_modes 1 --verbosity 2 --refine_step 5 --keep_nonpolar_H --local_only
</pre>pre>
<h2>Use DCU to carry out docking</h2>
<p>Prepare a docking scrript such as DC_batch_5.sh:</p>
<pre lang="shell">
unidock="/public/software/apps/unidock/install/bin/unidock"
export LD_LIBRARY_PATH=/public/software/apps/unidock/boost_1_72_0/build_sif/lib:$LD_LIBRARY_PATH
cd /work/home/achktgbwrc/cu2_unidock_evaluation
$unidock --ligand_index /work/home/achktgbwrc/index/DC_batch_5.index --config dock.conf --dir DC_batch_5 --search_mode fast
</pre>
<p>Prepare a SLURM script such as DC_batch_5.slurm:</p>
<pre lang="shell">
#!/bin/bash
#SBATCH -J DC_batch_5
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p xahdtest
#SBATCH --gres=dcu:1

module purge
module load singularity/3.7.3

export BOOST_ROOT=/public/software/apps/unidock/boost_1_72_0/build_sif
export BOOST_INCLUDEDIR=$BOOST_ROOT/include
export BOOST_LIBRARYDIR=$BOOST_ROOT/lib
export SINGULARITYENV_LD_LIBRARY_PATH="${BOOST_ROOT}:\$LD_LIBRARY_PATH"

singularity exec -B /public:/public -B /work:/work /public/software/apps/DeepLearning/singularity/centos7.6-mpi4.0-gcc9.3-cmake3.21-make4.2-glibc2.29-glibcxx3.4.26-py3.8-dtk23.10.sif bash -c "source /opt/dtk/env.sh && source /opt/dtk/cuda/env.sh && bash DC_batch_5.sh"
</pre>
<p>Submit jobs:</p>
<pre lang="shell">
sbatch DC_batch_5.slurm
</pre>

