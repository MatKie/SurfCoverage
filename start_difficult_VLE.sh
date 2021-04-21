#!/bin/bash

temp=298.15
surfactants=$1
tot_surfactants=$(($1+$1))
name='SDS'
c_ion='NA'
dir='VLE_'${surfactants}

if [ -d $dir ];
then
    rm -r $dir
fi

cp -r template_VLE $dir
cd $dir
echo '*' > .gitignore

sed -i s/NRSURFACTANT/$surfactants/g topo/topol.top
echo ${name}' '${tot_surfactants} >> topo/topol.top

cd config 
sed -i s/NRSURFACTANT/$surfactants/g water_SDS.inp
sed -i s/NRSURFACTANT/$surfactants/g water_SDS_difficult.inp
/home/mk8118/packmol/packmol < water_SDS_difficult.inp >> packmol.out


gmx editconf -f water_SDS.pdb -box 8 8 40 -o water_SDS.gro
gmx grompp -f dummy.mdp -c water_SDS.gro -p ../topo/topol.top 
gmx genion -s topol.tpr -p ../topo/topol.top -o water_SDS_neutral.gro -neutral -pname $c_ion << EOF
SOL
EOF

gmx grompp -f dummy.mdp -c water_SDS_neutral.gro -p ../topo/topol.top 

gmx make_ndx -f topol.tpr << EOF
t W2
t SO4V9
t CM
t CT
t NA+
t CM | t CT
q
EOF

cd ../em
gmx grompp -f em_difficult.mdp -c ../config/water_SDS_neutral.gro -p ../topo/topol.top -n ../config/index.ndx
/usr/local/gromacs2019/gromacs2019_dp/bin/mdrun_d -v -s topol.tpr

cd ../eq_npt

gmx grompp -f npt_difficult.mdp -c ../em/confout.gro -p ../topo/topol.top -n ../config/index.ndx
mv GROMACS.sh SC_VLE_${surfactants}.sh

gmx mdrun -v -nsteps 10000 -s topol.tpr

echo dir >> started.list
