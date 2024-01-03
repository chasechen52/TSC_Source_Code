package models;

import io.jenetics.*;
import io.jenetics.engine.Engine;
import io.jenetics.engine.EvolutionResult;
import io.jenetics.engine.EvolutionStatistics;
import io.jenetics.util.Factory;
import tool.RandomGraphGenerator;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static io.jenetics.engine.EvolutionResult.toBestPhenotype;
import static io.jenetics.engine.Limits.bySteadyFitness;

public class GAModel {

    private int mServersNumber;   // 服务器数量
    private int[][] mAccessMatrix;   // 可访问矩阵,初始为全0矩阵
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
        @SuppressWarnings("unchecked")
        List<T> dest = (List<T>) in.readObject();
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
            }
            mDegrees.put(i, access); // access为广义上的度
        }
        System.out.println("各节点度数：" + mDegrees);
    }

    public static int getMapMaxValueKey(Map<Integer, Integer> map)  // 取出最大value对应的Key值,Int
    {
        List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer, Integer>>(map.entrySet());
        Collections.sort(list, (o1, o2) -> (o1.getValue().intValue() - o2.getValue().intValue()));
        int key = list.get(list.size() - 1).getKey();
        return key;
    }


    public int getMapMinValue(Map<Integer, Integer> maps)  // 取出最小value
    {
        Comparator<Map.Entry<Integer, Integer>> valCmp = new Comparator<Map.Entry<Integer, Integer>>() {
            @Override
            public int compare(Map.Entry<Integer, Integer> o1, Map.Entry<Integer, Integer> o2) {
                return o1.getValue().intValue() - o2.getValue().intValue();  // 升序排序
            }
        };
        List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer, Integer>>(maps.entrySet()); // 传入maps实体
        Collections.sort(list, valCmp);

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

    private static int[][] GraphGenerate(int serversNumber, double density) {
        RandomGraphGenerator graphGenerator = new RandomGraphGenerator(serversNumber, density);
        graphGenerator.createRandomGraph();
        int[][] adjacencyMatrix = graphGenerator.getRandomGraphAdjacencyMatrix();
        int[][] distanceMatrics = graphGenerator.getRandomGraphDistanceMatrix();

        // int [][] distanceMatrics = new int [][]
        // {
        // 	{0,2,2,2,2,2,2,2,1,1},
        // 	{2,0,2,2,2,2,2,2,1,1},
        // 	{2,2,0,2,2,2,2,2,1,1},
        // 	{2,2,2,0,2,2,2,2,1,1},
        // 	{2,2,2,2,0,2,2,2,1,1},
        // 	{2,2,2,2,2,0,2,2,1,1},
        // 	{2,2,2,2,2,2,0,2,1,1},
        // 	{2,2,2,2,2,2,2,0,1,1},
        // 	{1,1,1,1,1,1,1,1,0,1},
        // 	{1,1,1,1,1,1,1,1,1,0},

        // };


        System.out.println("\n------- AdjMatrix -------"); // 打印邻接矩阵方便绘图
        for (int i = 0; i < serversNumber; ++i) {
            System.out.println(Arrays.toString(adjacencyMatrix[i]) + ",");
        }
        System.out.println("\n");

        return distanceMatrics;
    }

    // 轮盘赌选择
    private static List<Genotype<BitGene>> rouletteWheelSelection(List<Genotype<BitGene>> population, Map<Integer, Integer> mDegrees) {
        int chromeosomeTotalDgree = mDegrees.values().stream().mapToInt(Integer::intValue).sum();
        double[] cumulativeProbabilities = new double[population.size()];

        // 计算累积概率
        double cumulativeProbability = 0;
        for (int i = 0; i < population.size(); i++) {
            cumulativeProbability += (double) mDegrees.get(i) / chromeosomeTotalDgree;
            cumulativeProbabilities[i] = cumulativeProbability;
        }

        // 进行轮盘赌选择
        return Stream.generate(() -> {
                    double randomValue = Math.random();
                    for (int i = 0; i < cumulativeProbabilities.length; i++) {
                        if (randomValue <= cumulativeProbabilities[i]) {
                            return population.get(i);
                        }
                    }
                    return null; // 这里不会执行到，只是为了符合 Stream.generate 的要求
                })
                .limit(population.size())
                .collect(Collectors.toList());
    }
    // Define the fitness function
    private static Integer fitness(Genotype<BitGene> genotype) {
        // Your fitness calculation logic goes here
        // For each server, calculate its degree (connectivity) and multiply it by the gene value
        // Sum up the results for all servers to get the fitness value
        int fitness = 0;
        for (int i = 0; i < 10; i++) {
            int degree = mDegrees.get(i);
            int geneValue = genotype.getChromosome().getGene(i).getBit() ? 1 : 0;
            fitness += degree * geneValue;
        }
        return fitness;
    }

    public void runGACost(int population, int serverNumber) {
        // 首先进行矩阵转化，计算度
        ConvertDistoAccAndCalDegree();
        // 1.) Define the genotype (factory) suitable for the problem.
        // 我们创建了长度为 s 的BitChromosome，染色体中包含 1 的概率等于 0.5。
        Factory<Genotype<BitGene>> gtf = Genotype.of(BitChromosome.of(serverNumber, 0.5));
        System.out.println("gtf" + gtf);
        Engine<BitGene, Integer> engine = Engine.builder(GAModel::fitness, gtf)
                .populationSize(population)
                .selector(new RouletteWheelSelector<>())
                .alterers(new Mutator<>(0.2))
                .build();
        EvolutionStatistics<Integer, ?> statistics = EvolutionStatistics.ofNumber();

        // Run the genetic algorithm
        Phenotype<BitGene, Integer> result = engine.stream()
                .limit(bySteadyFitness(7))
                .peek(statistics)
                .collect(EvolutionResult.toBestPhenotype());
        // Print the best solution
        System.out.println("statistics \n" + statistics);
        System.out.println("Best Solution: " + result);

    }
}
