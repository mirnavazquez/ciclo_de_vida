for i in $(ls -d */)
do
cd $i
#cp t_data.ctab t_data_original.ctab
#cp t_data.ctab  y.ctab
rm t_data.ctab 
grep [t,r]RNA y.ctab | cut -f 6 > rna_names.tab
grep [t,r]RNA y.ctab | cut -f 1-9 > y_1.ctab
grep [t,r]RNA y.ctab | cut -f 11- > y_2.ctab
paste y_1.ctab rna_names.tab y_2.ctab > new_rna.tab
grep -v [t,r]RNA y.ctab > no_rna.tab
cat no_rna.tab new_rna.tab | sort -n -k 1,1 > t_data.ctab
cd ..
done


