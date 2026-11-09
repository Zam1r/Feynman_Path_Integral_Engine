$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path build | Out-Null
g++ -O3 -std=c++17 -pthread -I cpp -o build/feynman.exe cpp/main.cpp
g++ -O3 -std=c++17 -mwindows -I cpp -I third_party/raylib/src cpp/studio.cpp third_party/raylib/src/libraylib.a -lopengl32 -lgdi32 -lwinmm -static-libgcc -static-libstdc++ "-Wl,-Bstatic" -lwinpthread "-Wl,-Bdynamic" -o build/PathIntegral.exe
$gpp = (Get-Command g++).Source
Copy-Item (Join-Path (Split-Path $gpp) "libwinpthread-1.dll") -Destination "build\libwinpthread-1.dll" -Force
if ($args -contains "suite") {
    .\build\feynman.exe suite data
    python python/plot_all.py
}
