本程序自动生成top文件

*力场文件
*gaussview生成GJF文件，目前只能识别xyz坐标系。
*手动添加match文件，给每一个原子赋予原子的力场名称。
*安照需要修改template文件，控制top文件的输出行为。

在control.txt中进行初始化操作。

双击main.exe,输出结果。


MATCH FILE


详细说明文档：
程序由三个模块组成，
第一步读取gjf里的结构信息，生成键长，键角，二面角等。
structure.c
文件为readgjf.c graphic.c 

第二部分读取力场数据库，生成基本数据视图。
conndb.c ffpara.c
readff.c ffdb.c 
ffdb2.c readdb2.c

数据视图操作函数在opffdb.c中。
第三部分将前两部分链接起来，输出itp文件。
communicate .c
handshake.c
输出位于print.c

另外还包含一些工具函数，位于tools和utils文件夹。


注意涉及到的文件
gjf文件 由gaussview生成，直角坐标系格式。必须
gjf_match.txt文件，为gaussiview生成的文件赋予力场原子名称，需手动添加。必须
ff文件 数据库文件，可根据需求，自行添加数据。必须
目前支持A & B两种数据ff
itp的template文件，控制itp最终输出的一些行为。必须
control.txt控制输入输出。必须

注意原子变化在高斯中要按charge group组织