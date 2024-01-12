package models;

import io.jenetics.*;
import io.jenetics.engine.*;
import io.jenetics.stat.DoubleMomentStatistics;

import java.io.*;
import java.util.*;

public class GAModel {

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

    private int mPopulation; // 种群数量

    public GAModel(int serversNumber, int[][] distancematrix, int hops) // 构造函数
    {
        mServersNumber = serversNumber;
        mDistanceMatrix = distancematrix;
        mhops = hops;

        mAccessMatrix = new int[mServersNumber][mServersNumber];
        mSelectedServerList = new ArrayList<>();
        mDegrees = new HashMap<>();
        mDataPacketsNeed = new HashMap<>();
    }


    public static <T> List<T> deepCopy(List<T> src) throws IOException, ClassNotFoundException {
        ByteArrayOutputStream byteOut = new ByteArrayOutputStream();
        ObjectOutputStream out = new ObjectOutputStream(byteOut);
        out.writeObject(src);

        ByteArrayInputStream byteIn = new ByteArrayInputStream(byteOut.toByteArray());
        ObjectInputStream in = new ObjectInputStream(byteIn);
        @SuppressWarnings("unchecked") List<T> dest = (List<T>) in.readObject();
        return dest;
    }

    public void initDataPacketsNeed(int requiredpackets) {
        for (int key = 0; key < mServersNumber; ++key) // 将mDataPacketsNeed重置
        {
            mDataPacketsNeed.put(key, requiredpackets);
        }
    }

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
                System.out.print(mAccessMatrix[i][j] + ",");
            }
            System.out.println();
            mDegrees.put(i, access); // access为广义上的度
        }
        // System.out.println("各节点度数：" + mDegrees);
    }

    public static int getMapMaxValueKey(Map<Integer, Integer> map)  // 取出最大value对应的Key值,Int
    {
        List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer, Integer>>(map.entrySet());
        list.sort(Comparator.comparingInt(Map.Entry::getValue));
        return list.get(list.size() - 1).getKey();
    }

    public int getMapMaxValue(Map<Integer, Integer> maps)  // 取出最小value
    {
        Comparator<Map.Entry<Integer, Integer>> valCmp = new Comparator<Map.Entry<Integer, Integer>>() {
            @Override
            public int compare(Map.Entry<Integer, Integer> o1, Map.Entry<Integer, Integer> o2) {
                return o2.getValue() - o1.getValue();  // 升序排序
            }
        };
        List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer, Integer>>(maps.entrySet()); // 传入maps实体
        list.sort(valCmp);

        return list.get(0).getValue();
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

    public int selectServerwithMostDegrees() // 选择有最高度的节点
    {
        return getMapMaxValueKey(mDegrees);
    }

    public boolean checkPacketsRequired() // 检查是否所有点都已满足，若是返回true
    {
        int count = 0;
        for (int val : mDataPacketsNeed.values()) {
            if (val == 0) {
                count++;
            }
        }
        if (count == mDataPacketsNeed.size()) {
            return true;
        }
        return false;
    }

    // Define the fitness function
    private static double fitness(Genotype<IntegerGene> PLSolution) {
        Chromosome<IntegerGene> M_number_chromosome = PLSolution.getChromosome(0);
        Chromosome<IntegerGene> data_placement_chromosome = PLSolution.getChromosome(1);
        double N = 0;
        int M = M_number_chromosome.getGene(0).intValue();
        double fitness;
        boolean isFeasibleSolution = true;
        int dataCost;
        // int dataBenefit = 1;
        int seversWithoutEnoughData = 0;

        // 校验该数据放置策略是否可行：判断每个服务器是否能访问到足够多的数据块
        // 计算N：总数据块数; 计算数据增益：该数据放置在该位置能服务的服务器数量
        for (int i = 0; i < mServersNumber; i++) {
            N += data_placement_chromosome.getGene(i).intValue();
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
                    boolean hasData = data_placement_chromosome.getGene(j).intValue() == 1;
                    if (canAccess && hasData) --mRequired;
                    if (mRequired == 0) break;
                }
                if (mRequired > 0) {
                    isFeasibleSolution = false;
                    ++seversWithoutEnoughData;
                    break;

                }
                // System.out.println(isFeasibleSolution);
                // System.out.println(mRequired);
            }
        }
        // dataCost = -N;
        // 适应度：å数据增益 + ß成本 + ∂是否可行解
        double a = 1;
        // double b = 2.5;
        // double c = 1;

        double b = 500;
        // double isFeasibleSolutionValue = isFeasibleSolution ? 1 : -1;
        // dataBenefit = 1 / dataBenefit;
        // double adaptationItem = a * dataCost + b * N + c *  1 / M;
        // 适应项：
        double adaptationItem = a * N * (1 + (double) 1 / M);
        // 惩罚项：不可行的节点数 / 总服务器数
        double penaltyItem = b * seversWithoutEnoughData / mServersNumber;
        fitness = adaptationItem + penaltyItem;
        if (isFeasibleSolution) {
            System.out.println("可行解： " + "M=" +M + ",  数据分布：" + data_placement_chromosome);
        }
        return -fitness;
    }
    // 为1的度取倒数 1/k ，度的倒数之和相加
    // 整个演化过程
    // 最优解和其他解的区别

    public void runGACost(int population, int serverNumber) {
        // 首先进行矩阵转化，计算度
        ConvertDistoAccAndCalDegree();

        IntegerChromosome data_placement_chromosome = IntegerChromosome.of(0, 1, serverNumber);
        IntegerChromosome M_number_chromosome = IntegerChromosome.of(2, getMapMinValue(mDegrees), 1);

        Genotype<IntegerGene> PLSolution = Genotype.of(M_number_chromosome, data_placement_chromosome);

        System.out.println("PLSolution" + PLSolution);

        Engine<IntegerGene, Double> engine = Engine.builder(GAModel::fitness, PLSolution)
                // .constraint(new SolutionConstraint())
                .populationSize(population)
                .offspringFraction(0.8)
                .survivorsFraction(0.2)
                .selector(new TournamentSelector<>()) // 选择
                .alterers(new Mutator<>(0.2))  // 变异概率0.2
                .build();

        EvolutionStatistics<Double, ?> Statistics = EvolutionStatistics.ofNumber();

        EvolutionResult<IntegerGene, Double> result = engine.stream()
                .limit(Limits.bySteadyFitness(100))
                .peek(Statistics)
                .collect(EvolutionResult.toBestEvolutionResult());

        double M =  result.getBestPhenotype().getGenotype().getChromosome(0).getGene(0).intValue();
        double N = 0;
        for (int i = 0; i < data_placement_chromosome.length(); i++) {
            N += result.getBestPhenotype().getGenotype().getChromosome(1).getGene(i).intValue();
        }

        // System.out.println("statistics \n" + Statistics);
        System.out.println("Best Solution: " + result.getBestPhenotype());
        System.out.println("GAModel Cost: " + N/M);

    }
}
