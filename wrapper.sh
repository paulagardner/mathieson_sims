#!/bin/bash
source activate mathieson_sims
python model.py

bash downstream.sh

python plotting.py