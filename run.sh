g++ src/Main.cpp -o main -Ofast -Wall -Wextra
./main
source env/bin/activate
python scripts/plotting2D.py output.txt