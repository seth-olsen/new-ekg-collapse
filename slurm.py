# for submitting bash scripts

from subprocess import call
s = " "

program_name = "ekg"
outfile = "ekg"
ic_Amp = "0.004"
resn0 = 16

script_name = "mr.slurm"
job_name = "r_04"
num_nodes = "1"
tasks_per_node = "1"
walltime = "24:00:00"
mem_lim = "16000MB"

cmd_line_args = "-write_xp 0" + s + \
                "-write_abp 0" + s + \
                "-write_mtot 0" + s + \
                "-write_maspect 0" + s + \
                "-write_outnull 0" + s + \
                "-write_ricci 0" + s \
                + \
                "-ic_Dsq 5.0" + s + "-ic_r0 7.5" + s \
                + \
                "-nresn 1" + s + \
                "-resn1 2" + s + "-resn2 1" + s \
                + \
                "-lastpt 500" + s + "-save_pt 1" + s + \
                "-nsteps 2000" + s + "-save_step 4" + s + \
                "-lam 0.25" + s \
                + \             
                "-rmin 0.0" + s + "-rmax 50.0" + s \
                + \
                "-dspn 0.5" + s + \
                "-tol 0.000000001" + s + \
                "-ell_tol 0.000000000001" + s + \
                "-ell_up_weight 0.5" + s + \
                "-maxit 100" + s \
                + \
                "-psi_hyp 1" + s + \
                "-horizon_search 0" + s \
                + \
                "-check_step 10" + s + \
                "-write_res 0" + s + \
                "-write_itn 0" + s + \
                "-write_ires_xp 0" + s + \
                "-write_ires_abp 0" + s + \
                "-dspn_bound 0" + s + \
                "-somm_cond 1" + s + \
                "-dspn_psi 0" + s + \
                "-dr3_up 0" + s + \
                "-clean_hyp 0" + s + \
                "-clean_ell 0" + s + \
                "-static_metric 0" + s + \
                "-zero_pi 0" + s + \
                "-resnorm_type 0" + s + \
                "-r2m 0.0" + s + \
                "-hold_const lambda" + s + \
                "-outfile " + outfile + s + \
                "-ic_Amp " + ic_Amp + s + \
                "-resn0 " + resn0

script_text = "#!/bin/sh" + "\n\n" + \
              "#SBATCH -N " + num_nodes + "\n" + \
              "#SBATCH --ntasks-per-node=" + tasks_per_node + "\n" + \
              "#SBATCH -J " + job_name + "\n" + \
              "#SBATCH --time=" + walltime + "\n" + \
              "#SBATCH --mem=" + mem_lim + "\n" + \
              "#SBATCH --mail-user=srolsen@princeton.edu" + "\n" + \
              "#SBATCH --mail-type=begin" + "\n" + \
              "#SBATCH --mail-type=end" + "\n\n" + \
              "./" + program_name + s + cmd_line_args + "\n"

f = open(script_name, 'w')
f.write(script_text)
f.close()

bash_text = "sbatch ./" + script_name
call(bash_text, shell=True)


