./compile.sh
./build/lbm data/ball.txt 5 -gpu
source env/bin/activate
python scripts/plotting2D.py output.txt