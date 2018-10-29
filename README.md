# Dsp_IIR
transform of from IIR analog filter  to IIR digital filter

- ## 环境：
  - python-2.7
  - mpmath-1.0.0
  - sympy-1.3.0
- ## 环境配置方法：
  - win：打开命令行，cd 进入python安装目录下的Scripts，将项目中的两个whl文件复制进去，输入的命令在“捕获”中，如没有pip先使用easy安装pip
  - Ubuntu：pip install即可
- ## 主要功能：
  - 模拟低通到模拟高通、带通、带阻的转换
  - Analog Low Pass Filter to Analog High Pass/Band Pass/Bandstop Filter
  - 模拟低通到数字高通、带通、带阻的转换，直接法
  - Analog Low Pass Filter to Digital High Pass/Band Pass/Bandstop Filter(Way I)
  - 模拟低通到模拟高通、带通、带阻的转换，双线性法
  - Analog Low Pass Filter to Digital High Pass/Band Pass/Bandstop Filter(Way II)
  - 可以点阵绘图，格式化输出
- ## 主要参数：
  - isChebyOrButt：1表示切比雪夫，0表示巴特沃斯
  - f_p, f_st 分别表示通带阻带
  - 输入函数的单位需要注意，提供了两个转换函数
    - 频率->数字频率
    - 频率->角频率
- ## 声明：
 本文件是个人短时间之内制作，难免有很多的BUG，欢迎测试和Report
