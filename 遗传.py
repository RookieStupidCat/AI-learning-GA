import matplotlib.pyplot as plt 
import math 
import random 
 
""" 
本项目参考自https://blog.csdn.net/hiudawn/article/details/80144221的博客
在此基础上做了一些优化修改
本例求解的目标表达式为： 
y = x * sin(10 * pi * x) + 1.0 , x->[-1,2] 

根据情况改变plot_obj_func(x)里的表达式
以及downer_limit，upper_limit所代表的自变量范围
""" 
# 目标函数 
def plot_obj_func(x): 
    return x * math.sin(10 * math.pi * x) + 1  

 
def main(): 
    pop_size = 500  # 种群数量 
    upper_limit = 2  # 基因中允许出现的最大值 ————这里根据自变量范围修改
    downer_limit = -1   #自变量下限
    chromosome_length = 10  # 染色体长度 ————跟精读和概率有关
    iter = 500 # 遗传代数
    pc = 0.6 # 杂交概率 
    pm = 0.01  # 变异概率 
    results_X = []  # 存储每一代的最优解的基因序列
    results_Y = []  # 存储每一代的最优解的个体适应度
    pop = init_population(pop_size, chromosome_length) # 随机初始化染色体，记录在pop[[],[],[]...]中
    for i in range(iter): 
        # 个体评价，本函数中就是f(x)
        obj_value = calc_obj_value(pop, chromosome_length, upper_limit,downer_limit)   
        # 个体适应度，达不到筛选条件的归0 
        fit_value = calc_fit_value(obj_value)  
        # 第一个是最优基因序列, 第二个是对应的最佳个体适度 
        best_individual, best_fit = find_best(pop, fit_value)  
        # 下面这句就是存放每次迭代的最优x值是最佳y值 
        results_X.append(binary2decimal(best_individual, upper_limit,downer_limit, chromosome_length)) 
        results_Y.append(best_fit)

        # 选择————轮盘赌 
        selection(pop, fit_value)
        # 染色体交叉（最优个体之间进行0、1互换）  
        crossover(pop, pc , fit_value)   
        # 染色体变异（其实就是随机一个进行0、1取反） 
        mutation(pop, pm , fit_value)  

    print('y = x * sin(10 * pi * x) + 1.0') # 改动函数此处需更改
    print("x = %f, y = %f" % (results_X[-1], results_Y[-1])) 
    # 看种群点的选择 
    plt.scatter(results_X,results_Y, s=3, c='r') 
    X1 = [i / float(100) for i in range(-100, 200, 1)] 
    Y1 = [plot_obj_func(x) for x in X1] 
    plt.plot(X1, Y1) 
    plt.show() 
    # 查看迭代曲线 
    plot_iter_curve(iter,results_Y) 
 
 
 
# 看看当前种群个体的落点情况 
def plot_currnt_individual(X, Y): 
    X1 = [i / float(100) for i in range(-100, 200, 1)] 
    Y1 = [plot_obj_func(x) for x in X1] 
    plt.plot(X1, Y1) 
    plt.scatter(X, Y, c='r', s=5) 
    plt.show() 
 
 
# 看看最终的迭代变化曲线 
def plot_iter_curve(iter,results_Y): 
    X = [i for i in range(iter)] # 迭代次数为横坐标
    # Y = [results_X[i] for i in range(iter)] # 每次迭代的最优值
    plt.plot(X, results_Y) 
    plt.show() 
 
 
# 计算2进制序列代表的数值 
def binary2decimal(binary, upper_limit,downer_limit, chromosome_length): 
    t = 0 
    for j in range(len(binary)): 
        t += binary[j] * 2 ** j 
    t = downer_limit + t * (upper_limit - downer_limit) / (2 ** chromosome_length - 1) 
    return t 
 
# 初始化基因（二进制 
def init_population(pop_size, chromosome_length): 
    # 形如[[0,1,..0,1],[0,1,..0,1]...] 
    pop = [[random.randint(0, 1) for i in range(chromosome_length)] for j in range(pop_size)] 
    return pop 
 
 
# 解码并计算值 
def decode_chromosome(pop, chromosome_length, upper_limit,downer_limit): 
    X = [] 
    for ele in pop: 
        temp = 0 
        # 二进制变成实数，种群中的每个个体对应一个数字 
        for i, coff in enumerate(ele): 
            # 就是把二进制转化为十进制的 
            temp += coff * (2 ** i) 
 
        # 这个是把前面得到的那个十进制的数，再次缩放为另一个实数 
        # 注意这个实数范围更广泛，可以是小数了，而前面二进制解码后只能是十进制的数 
        # 参考https://blog.csdn.net/robert_chen1988/article/details/79159244 
        X.append(downer_limit + temp * (upper_limit - downer_limit) / (2 ** chromosome_length - 1)) 
        # print(ele)
    return X 
 
#  计算种群个体优劣
def calc_obj_value(pop, chromosome_length, upper_limit,downer_limit): 
    obj_value = [] 
    # 先解码成十进制
    X = decode_chromosome(pop, chromosome_length, upper_limit,downer_limit) 
    for x in X: 
        # 把缩放过后的那个数，带入我们要求的公式中 
        # 种群中个体有几个，就有几个这种“缩放过后的数” 
        obj_value.append(plot_obj_func(x)) 
    # 这里先返回带入公式计算后的数值列表，作为种群个体优劣的评价 
    return obj_value 
 
 
# 淘汰 
def calc_fit_value(obj_value): 
    fit_value = [] 
    # 更改c_min会改变淘汰的下限 
    # 但是如果设置过大，有可能影响了全局最优的搜索 
    # 本例用最简单的常数筛选，更改后需据情况修改
    c_min = 2 
    for value in obj_value: 
        if value > c_min: 
            temp = value 
        else: 
            temp = 0. 
        fit_value.append(temp) 
    # fit_value保存的是活下来的值 
    return fit_value 
 
 
# 找出最优解和最优解的基因编码 
def find_best(pop, fit_value): 
    # 用来存最优基因编码 
    best_individual = [] 
    # 先假设第一个基因的适应度最好 
    best_fit = fit_value[0] 
    for i in range(1, len(pop)): 
        if (fit_value[i] > best_fit): 
            best_fit = fit_value[i] 
            best_individual = pop[i] 
    # best_fit是值 
    # best_individual是基因序列 
    return best_individual, best_fit 
 
 
# 计算概率和 
def cum_sum(fit_value): 
    # 输入[1, 2, 3, 4, 5]，返回[1,3,6,10,15]
    # 这个地方遇坑，局部变量如果赋值给引用变量，在函数周期结束后，引用变量也将失去这个值 
    temp = fit_value[:] 
    for i in range(len(temp)): 
        fit_value[i] = (sum(temp[:i + 1])) 
 
 
# 轮盘赌选择 
def selection(pop, fit_value): 
    # 存储转化后的概率值
    p_fit_value = [] 
    # 适应度总和 
    total_fit = sum(fit_value) 
    # 归一化，转化为概率 
    for i in range(len(fit_value)): 
        p_fit_value.append(fit_value[i] / total_fit) 
    # 概率和排序 
    cum_sum(p_fit_value) 
    pop_len = len(pop) 
    # 为每一个染色体转一次轮盘，记录结果到ms[]
    ms = sorted([random.random() for i in range(pop_len)]) 
    fitin = 0 #两个类似指针的作用
    newin = 0 
    newpop = pop[:] 
    # 转轮盘选择法 
    '''这里效率有点低，试试优化成二分查找'''
    while newin < pop_len: 
        # 找到与随机到的概率相符的染色体区间/序号
        if(p_fit_value[newin]!=p_fit_value[newin-1]):# 不处理0——以淘汰的个体
            if (ms[newin] < p_fit_value[fitin]): # 找到
                newpop[newin] = pop[fitin] #用随机到的染色体替换掉newin位置的染色体
                newin = newin + 1 
                fitin = 0
            else: 
                fitin = fitin + 1 
        else:
            newin = newin + 1
    # 这里注意一下，因为random.random()不会大于1，所以保证这里的newpop规格会和以前的一样 
    # 而且这个pop里面会有不少重复的个体，保证种群数量一样 
    pop = newpop[:] 
 
 
# 杂交 
def crossover(pop, pc,fit_value): 
    # 一定概率杂交，主要是杂交种群种相邻的两个个体 
    pop_len = len(pop) 
    for i in range(pop_len - 1): 
        # 随机看看达到杂交概率没 
        if fit_value[i]!=0:
            pass
        elif ( random.random() < pc): 
            # 随机选取杂交点，然后交换数组 
            cpoint = random.randint(0, len(pop[0])) 
            temp1 = [] 
            temp2 = [] 
            temp1.extend(pop[i][0:cpoint]) 
            temp1.extend(pop[i + 1][cpoint:len(pop[i])]) 
            temp2.extend(pop[i + 1][0:cpoint]) 
            temp2.extend(pop[i][cpoint:len(pop[i])]) 
            pop[i] = temp1[:] 
            pop[i + 1] = temp2[:] 
 
 
# 基因突变 
def mutation(pop, pm,fit_value): 
    px = len(pop) 
    py = len(pop[0]) 
    # 每条染色体随便选一个杂交 
    for i in range(px): 
        if fit_value[i]!=0:
            pass
        elif  (random.random() < pm): 
            mpoint = random.randint(0, py - 1) 
            if (pop[i][mpoint] == 1): 
                pop[i][mpoint] = 0 
            else: 
                pop[i][mpoint] = 1 
 
 
if __name__ == '__main__': 
    main() 
 