#!/bin/bash
source activate ofup
export PATH=/code/:$PATH
jupyter notebook --ip 0.0.0.0 --no-browser --allow-root
