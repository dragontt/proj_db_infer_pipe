import numpy

fn_chip = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_physical_network_chip/chip_net_genome_wide.adjmtr'
fn_pwm = '/home/mblab/ykang/proj_db_infer_pipe/resources/fly_net/resources/networks/physical/motif_net.adjmtr'

# get subnet based on pwm nonzero gids
chip = numpy.loadtxt(fn_chip)
pwm = numpy.loadtxt(fn_pwm)
# gid_indices = numpy.intersect1d(numpy.unique(numpy.nonzero(chip)[1]), numpy.unique(numpy.nonzero(pwm)[1]))
gid_indices = numpy.unique(numpy.nonzero(pwm)[1])
chip_sub = chip[:,gid_indices]
pwm_sub = pwm[:,gid_indices]

# overlap of rids
rid_indices = numpy.intersect1d(numpy.unique(numpy.nonzero(chip_sub)[0]), numpy.unique(numpy.nonzero(pwm_sub)[0]))
chip_sub = chip_sub[rid_indices,:]
pwm_sub = pwm_sub[rid_indices,:]
overlap_sub = numpy.multiply(chip_sub, pwm_sub)

# compute stats
total_chip_sub_positive = len(numpy.nonzero(chip_sub)[0])
total_pwm_sub_positive = len(numpy.nonzero(pwm_sub)[0])
total_overlap_interaction = len(numpy.nonzero(overlap_sub)[0])
size_pwm_sub = pwm_sub.shape[0]*pwm_sub.shape[1]

print "ChIP on PWM: ", (float(total_chip_sub_positive)/969, float(total_overlap_interaction)/size_pwm_sub)
print "chance on PWM : ", float(total_pwm_sub_positive)/size_pwm_sub
