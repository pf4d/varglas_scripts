import sys
from varglas.data.data_factory    import DataFactory
from varglas.utilities            import DataInput, DataOutput

thklim = 200.0

# collect the raw data :
searise  = DataFactory.get_searise(thklim = thklim)
rignot   = DataFactory.get_gre_rignot_updated()

# create data objects to use with varglas :
dsr     = DataInput(searise)
drg     = DataInput(rignot)

# change the projection of the measures data to fit with other data :
drg.change_projection(dsr)

# inspect the data values :
do    = DataOutput('results_pre/')
do.write_one_file('vmag',           drg.get_projection('U_ob'))
#do.write_one_file('h',              dbm.get_projection('H'))
#do.write_one_file('Ubmag_measures', dbv.get_projection('Ubmag_measures'))
#do.write_one_file('sr_qgeo',        dsr.get_projection('q_geo'))
#exit(0)


