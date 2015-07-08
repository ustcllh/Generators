NAME=$1
g++ $NAME -Werror -Wall -O2 -o "${NAME/%.C/}.exe" `/Users/cfmcginn/Projects/Generators/PYTHIA/pythia-install/bin/pythia8-config --cxxflags --libs` -pthread -std=c++11 -m64 -I/opt/local/libexec/root5/include/root -L/opt/local/libexec/root5/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread  -lm -ldl

#big string after `` is root-config --cflags --libs minus the --stdlib part, for some reason it fucking hates that