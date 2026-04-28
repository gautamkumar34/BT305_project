#!/bin/bash
# run_train.sh

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES

echo "Starting hERG Chemprop Training Pipeline..."
PYTHONPATH=. venv/bin/python3 -c "
from src.tox_model import train_and_save_model
train_and_save_model('data/herg_training_final.csv', 'models/chemprop_herg', epochs=30)
"
