1.在目录project下输入 make 会新建一个目录bin，会生成可执行文件到目录bin中。
2.在目录project下输入 make run 会运行make所得到的可执行文件（该指令并不会编译，所以一定要先make）并将结果绘图。
3.make clean 可以将运行结果（result目录中的所有文件以及fig目录中的.jpg文件）以及可执行文件清除。
4.样条的实现在./src/include/splines.h头文件中
5.参数以json文件格式输入，具体见parameter.json，内有注释（@开头的即注释）。
6.测试结果是report.pdf
7.设计文档分为设计思路（ 设计文档.pdf）和由Doxygen生成的代码文档（有pdf版本(design document.pdf)和网页版本，网页版本进入目录design document html打开annotated.html即可）
