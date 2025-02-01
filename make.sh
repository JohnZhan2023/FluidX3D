#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#partition name
#SBATCH --partition=viscam,svl
#################
#number of GPUs
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --account=viscam
#################
#set a job name
#SBATCH --job-name="neural_wind_tunnel"
#################
#a file for job output, you can check job progress, append the job ID with %j to make it unique
##SBATCH --output=test_output/%j.out
#################
# a file for errors from the job
##SBATCH --error=test_output/%j.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm, 144
#SBATCH --time=12:00:00
#################
# Quality of Service (QOS); think of it as sending your job into a special queue; --qos=long for with a max job length of 7 days.
# uncomment ##SBATCH --qos=long if you want your job to run longer than 48 hours, which is the default for normal partition,
# NOTE- in the hns partition the default max run time is 7 days , so you wont need to include qos, also change to normal partition
# since dev max run time is 2 hours.
##SBATCH --qos=long
# We are submitting to the dev partition, there are several on sherlock: normal, gpu, bigmem (jobs requiring >64Gigs RAM)
##SBATCH -p dev
#################
# --mem is memory per node; default is 4000 MB per CPU, remember to ask for enough mem to match your CPU request, since
# sherlock automatically allocates 4 Gigs of RAM/CPU, if you ask for 8 CPUs you will get 32 Gigs of RAM, so either
# leave --mem commented out or request >= to the RAM needed for your CPU request.  It will also accept mem. in units, ie "--mem=4G"
#SBATCH --mem=80G
#################
# Have SLURM send you an email when the job ends or fails, careful, the email could end up in your clutter folder
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email
#SBATCH --mail-user=jiahao@cs.stanford.edu
# list out some useful information
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NAME="$SLURM_JOB_NAME
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR
#now run normal bash commands


#!/usr/bin/env bash
# command line argument(s) for make.sh: device ID(s); if empty, FluidX3D will automatically choose the fastest available device(s)

case "$(uname -a)" in # automatically detect operating system and X11 support on Linux
	 Darwin*) target=macOS                                                       ;;
	*Android) target=Android                                                     ;;
	 Linux* ) if xhost >&/dev/null; then target=Linux-X11; else target=Linux; fi ;;
	*       ) target=Linux                                                       ;;
esac

#target=Linux-X11 # manually set to compile on Linux with X11 graphics
#target=Linux     # manually set to compile on Linux (without X11)
#target=macOS     # manually set to compile on macOS (without X11)
#target=Android   # manually set to compile on Android (without X11)

echo -e "\033[92mInfo\033[0m: Detected Operating System: "${target}
echo_and_execute() { echo "$@"; "$@"; }
if command -v make &>/dev/null; then # if make is available, compile FluidX3D with multiple CPU cores
	echo -e "\033[92mInfo\033[0m: Compiling with "$(nproc)" CPU cores."
	make ${target} -j$(nproc) # compile FluidX3D with makefile
else # else (make is not installed), compile FluidX3D with a single CPU core
	echo -e "\033[92mInfo\033[0m: Compiling with 1 CPU core. For faster multi-core compiling, install make with \"sudo apt install make\"."
	mkdir -p bin # create directory for executable
	rm -rf temp bin/FluidX3D # prevent execution of old executable if compiling fails
	case "${target}" in
		Linux-X11) echo_and_execute g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -O -Wno-comment -I./src/OpenCL/include -L./src/OpenCL/lib -lOpenCL -I./src/X11/include -L./src/X11/lib -lX11 -lXrandr ;;
		Linux    ) echo_and_execute g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -O -Wno-comment -I./src/OpenCL/include -L./src/OpenCL/lib -lOpenCL                                                    ;;
		macOS    ) echo_and_execute g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -O -Wno-comment -I./src/OpenCL/include -framework OpenCL                                                              ;;
		Android  ) echo_and_execute g++ src/*.cpp -o bin/FluidX3D -std=c++17 -pthread -O -Wno-comment -I./src/OpenCL/include -L/system/vendor/lib64 -lOpenCL                                                ;;
	esac
fi

if [[ $? == 0 ]]; then bin/FluidX3D "$@"; fi # run FluidX3D only if last compilation was successful
