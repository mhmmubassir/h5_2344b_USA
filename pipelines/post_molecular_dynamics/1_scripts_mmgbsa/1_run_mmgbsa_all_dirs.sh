#!/bin/bash
###############################################################################
#  run_all_MMGBSA.sh
#  ----------------------------------------------------------
#  • Walk through each first-level subdirectory
#  • Skip y.input_files/ and z.template/
#  • Create   <subdir>/mmgbsa/     if it doesn’t exist
#  • Copy CPLX.parm7 + produ.nc + template files
#  • Run cpptraj to split receptor / ligand topologies
#  • Submit the GBneck2 job with submit_MMGBSA.sh
###############################################################################

# Absolute path to the current working directory
current_dir=$(pwd)

# Template files that will be copied into every mmgbsa/ folder
files_to_copy=("cpptraj.in" "mmgbsa.in" "submit_MMGBSA.sh")

echo
echo "===================================================================="
echo " MMGBSA dispatcher started in: ${current_dir}"
echo " Will copy: ${files_to_copy[*]}"
echo "===================================================================="
echo

# ---------------------------------------------------------------------
# Loop through each immediate sub-directory
# ---------------------------------------------------------------------
for subdir in */ ; do

    # Skip template or input-file directories
    if [[ "${subdir}" == "y.input_files/" || "${subdir}" == "z.template/" ]]; then
        echo ">>> Skipping directory: ${subdir}"
        continue
    fi

    echo ">>> Checking directory: ${subdir}"
    full_subdir="${current_dir}/${subdir}"

    # Confirm that subdir really exists and is a directory
    if [[ -d "${full_subdir}" ]]; then
        echo "    ✔  Directory exists: ${full_subdir}"

        # Define workspace path
        workdir="${full_subdir}mmgbsa/"

        # Create workspace if it does not exist
        if [[ ! -d "${workdir}" ]]; then
            echo "    •  Creating workspace: ${workdir}"
            mkdir -p "${workdir}" || { echo "    ✖  mkdir failed"; exit 1; }
        else
            echo "    •  Workspace already present: ${workdir}"
        fi

        # ---------------- Copy main trajectory & topology ----------------
        cp -v "${full_subdir}CPLX.parm7"  "${workdir}" \
            && echo "    •  Copied CPLX.parm7"
        cp -v "${full_subdir}produ.nc"    "${workdir}" \
            && echo "    •  Copied produ.nc"

        # ---------------- Copy template input files ---------------------
        for file in "${files_to_copy[@]}"; do
            cp -v "${current_dir}/${file}" "${workdir}" \
                && echo "    •  Copied ${file}"
        done

        # ---------------- Enter workspace and launch jobs ---------------
        cd "${workdir}"
        if [[ $? -ne 0 ]]; then
            echo "    ✖  Cannot cd into ${workdir} — aborting"
            exit 1
        fi
        echo "    •  Entered workspace: $(pwd)"

        # Run cpptraj to generate receptor & ligand topologies
        cpptraj -p CPLX.parm7 -i cpptraj.in
        if [[ $? -eq 0 ]]; then
            echo "    •  cpptraj split completed successfully"
        else
            echo "    ✖  cpptraj error — skipping SLURM submission"
            cd "${current_dir}"; continue
        fi

        # Submit the MMGBSA SLURM job
        sbatch submit_MMGBSA.sh
        if [[ $? -eq 0 ]]; then
            echo "    •  SLURM job submitted"
        else
            echo "    ✖  sbatch submission failed"
        fi

        # Return to original directory
        cd "${current_dir}" || { echo "    ✖  Return cd failed"; exit 1; }
        echo
    else
        echo "    ✖  ${full_subdir} is not a directory — skipping."
        echo
    fi
done

echo "===================================================================="
echo " MMGBSA dispatcher finished."
echo "===================================================================="
