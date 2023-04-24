这是一个将分子动力学模拟和化学反应（Gillespie算法）耦合起来的包
示例如下


在打开julia之前，设置线程数，用来多线程加速，也可以不设置
通过在终端中输入 julia --threads auto 将线程数设置为本地CPU线程数(Julia1.7)；
也可以 julia --threads 4 指定为整数(Julia1.5)；
也可以设置环境变量，
Bash (Linux/macOS):
export JULIA_NUM_THREADS=4
C shell on Linux/macOS, CMD on Windows:
set JULIA_NUM_THREADS=4
Powershell on Windows:$env:
JULIA_NUM_THREADS=4


可以通过 Threads.nthreads() 查看当前线程数
然后将包导入
include("Mol_Gillespie.jl")

然后可以参考已经设置好的例子文件test.jl
include("test.jl")