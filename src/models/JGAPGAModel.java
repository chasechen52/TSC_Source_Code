package models;

import io.jenetics.IntegerChromosome;
import org.jgap.*;
import org.jgap.impl.DefaultConfiguration;
import org.jgap.impl.IntegerGene;

import java.util.*;


public class JGAPGAModel {

    private static int mServersNumber;   // 服务器数量
    private static int[][] mAccessMatrix;   // 可访问矩阵,初始为全0矩阵
    private int[][] mDistanceMatrix;  // 最短距离矩阵
    private int mhops;     // 允许的跳数
    private int mPacketsNeed;  // 最后需要的分块数

    private double mCost; // Cost
    private double mReplicaCost;
    private static Map<Integer, Integer> mDegrees; // 度
    private Map<Integer, Integer> mDataPacketsNeed; // 尚需数据块的个数
    private List<Integer> mSelectedServerList; // 数据节点列表

    private static ArrayList<List> mCplexSolutionList;

    private int mPopulation; // 种群数量


    public void ConvertDistoAccAndCalDegree() // 根据hops，将能访问到的节点全部置1，否则置0;同时计算度
    {
        for (int i = 0; i < mServersNumber; ++i) {
            int access = 0; // 可访问节点数
            for (int j = 0; j < mServersNumber; ++j) {
                if (mDistanceMatrix[i][j] <= mhops) // 当最短距离小于hops时，可访问节点数自增
                {
                    access++;
                    mAccessMatrix[i][j] = 1;
                } else mAccessMatrix[i][j] = 0;
                // System.out.print(mAccessMatrix[i][j] + ",");
            }
            // System.out.println();
            mDegrees.put(i, access); // access为广义上的度
        }
        // System.out.println("各节点度数：" + mDegrees);
    }

    public int getMapMinValue(Map<Integer, Integer> maps)  // 取出最小value
    {
        Comparator<Map.Entry<Integer, Integer>> valCmp = new Comparator<Map.Entry<Integer, Integer>>() {
            @Override
            public int compare(Map.Entry<Integer, Integer> o1, Map.Entry<Integer, Integer> o2) {
                return o1.getValue() - o2.getValue();  // 升序排序
            }
        };
        List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer, Integer>>(maps.entrySet()); // 传入maps实体
        list.sort(valCmp);

        return list.get(0).getValue();
    }

    public JGAPGAModel(int serversNumber, int[][] distancematrix, int hops) // 构造函数
    {
        mServersNumber = serversNumber;
        mDistanceMatrix = distancematrix;
        mhops = hops;

        mAccessMatrix = new int[mServersNumber][mServersNumber];
        mSelectedServerList = new ArrayList<>();
        mDegrees = new HashMap<>();
    }

    // 初始化数组的方法

    // 自定义适应度函数
    public static class MyFitnessFunction extends FitnessFunction {
        @Override
        protected double evaluate(IChromosome strategy_chromosome) {
            // System.out.println("chromosome: " + strategy_chromosome);
            // 在这里定义适应度函数的计算逻辑

            // 获取M的值
            int M = (int) strategy_chromosome.getGene(0).getAllele();
            double N = 0;
            double fitness;
            boolean isFeasibleSolution = true;
            int seversWithoutEnoughData = 0;

            // 校验该数据放置策略是否可行：判断每个服务器是否能访问到足够多的数据块
            // 计算N：总数据块数; 计算数据增益：该数据放置在该位置能服务的服务器数量
            for (int i = 0; i < mServersNumber; i++) {
                N += (int) strategy_chromosome.getGene(i).getAllele();
                // dataBenefit += data_placement_chromosome.getGene(i).intValue() * mDegrees.get(i);
            }
            if (N < M) {
                // 如果总数据块不够则直接判定不可行
                isFeasibleSolution = false;
            } else {
                // 总数据块数量足够，开始判断每个服务器需要的数据块是否足够
                for (int i = 0; i < mServersNumber; i++) {
                    int mRequired = M;
                    for (int j = 0; j < mServersNumber; j++) {
                        boolean canAccess = mAccessMatrix[i][j] == 1;
                        boolean hasData = (int) strategy_chromosome.getGene(j).getAllele() == 1;
                        if (canAccess && hasData) --mRequired;
                        if (mRequired == 0) break;
                    }
                    if (mRequired > 0) {
                        isFeasibleSolution = false;
                        ++seversWithoutEnoughData;
                        break;

                    }
                }
            }

            double a = 1;
            double b = 10000;
            // 适应项：
            double adaptationItem = a * N * (1 + (double) 1 / M);
            // 惩罚项：不可行的节点数 / 总服务器数
            double penaltyItem = b * seversWithoutEnoughData / mServersNumber;
            fitness = adaptationItem + penaltyItem;
            if (isFeasibleSolution) {
                System.out.println("可行解： " + "M=" +M + ",  数据分布：" + strategy_chromosome);
            }
            return 1 / fitness;
            // return 0;
        }
    }


    public void runGACost(int populationSize, int serverNumber, ArrayList<List> cplexSolutionList) throws InvalidConfigurationException {

        ConvertDistoAccAndCalDegree();
        Configuration configuration = new DefaultConfiguration();

        mCplexSolutionList = cplexSolutionList;

        IChromosome[] cplex_solution_chromosome_list = new IChromosome[cplexSolutionList.size()];

        for (int i = 0; i < cplex_solution_chromosome_list.length; i++) {
            List solution = cplexSolutionList.get(i);
            Gene[] genes = new Gene[serverNumber + 1];
            for (int j = 0; j < serverNumber + 1; j++) {
                if (j == 0) {
                    genes[j] = new IntegerGene(configuration, 1, getMapMinValue(mDegrees));
                    genes[j].setAllele(i+1);
                } else {
                    genes[j] = new IntegerGene(configuration, 0, 1);
                    genes[j].setAllele(0);
                }
            }
            for (Object index : solution) {
                int mIndex = (int) index;
                if (mIndex >= 0 && mIndex < serverNumber + 1) {
                    // genes[mIndex] = new IntegerGene(configuration, 0, 1);
                    genes[mIndex].setAllele(1);
                }
            }
            IChromosome feasible_solution = new Chromosome(configuration, genes);
            cplex_solution_chromosome_list[i] = feasible_solution;
        }

        // 创建一个遗传算法配置

        // 创建第一个染色体（长度为1，取值范围为0到10的整数）
        Gene[] genes = new Gene[serverNumber + 1];

        // 创建复合基因并将两个基因合并成一个染色体
        for (int i = 0; i < genes.length; i++) {
            if (i == 0) genes[i] = new IntegerGene(configuration, 1, getMapMinValue(mDegrees));
            else genes[i] = new IntegerGene(configuration, 0, 1);
        }

        // 创建染色体并将多个基因添加到其中
        IChromosome data_placement_strategy_chromosome = new Chromosome(configuration, genes);

        configuration.setSampleChromosome(data_placement_strategy_chromosome);
        configuration.setPopulationSize(populationSize);

        // 设置适应度函数
        FitnessFunction fitnessFunction = new MyFitnessFunction();
        configuration.setFitnessFunction(fitnessFunction);


        // 使用Genotype.randomInitialGenotype()创建一个随机初始化的个体的基因型
        Genotype genotype = Genotype.randomInitialGenotype(configuration);

        Population population = genotype.getPopulation();
        for (int number = 0; number < 10 ;number++) {
            for (int i = 0; i < cplex_solution_chromosome_list.length; i++) {
                population.addChromosome(cplex_solution_chromosome_list[i]);

            }
        }


        // 执行一定数量的进化操作
        int numberOfEvolutions = 100; // 你可以根据需要设置进化的次数
        for (int i = 0; i < numberOfEvolutions; i++) {
            genotype.evolve();
        }
        // 获取最优染色体
        IChromosome fittestChromosome = genotype.getFittestChromosome();

        double cost;
        double N = 0;

        for (int i = 1; i < 21; i++) {
            N += (int) fittestChromosome.getGene(i).getAllele();
        }
        int M = (int) fittestChromosome.getGene(0).getAllele();

        cost = N / M;

        // 打印最优染色体的基因值
        System.out.println("mGACost: " + cost);
    }
}




