package models;

import io.jenetics.*;
import io.jenetics.engine.*;
import sun.plugin2.os.windows.OSVERSIONINFOA;

import java.io.*;
import java.util.*;

public class GAModel {

    private int mServersNumber;   // 服务器数量
    private static int[][] mAccessMatrix;   // 可访问矩阵,初始为全0矩阵
    private int[][] mDistanceMatrix;  // 最短距离矩阵
    private int mhops;     // 允许的跳数
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
                // System.out.print(mAccessMatrix[i][j] + ",");
            }
            // System.out.println();
            mDegrees.put(i, access); // access为广义上的度
        }
        // System.out.println("各节点度数：" + mDegrees);
    }

    public static int getMapMaxValueKey(Map<Integer, Integer> map)  // 取出最大value对应的Key值,Int
    {
        List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer, Integer>>(map.entrySet());
        Collections.sort(list, (o1, o2) -> (o1.getValue().intValue() - o2.getValue().intValue()));
        int key = list.get(list.size() - 1).getKey();
        return key;
    }

    public static int getMapMaxValue(Map<Integer, Integer> maps)  // 取出最小value
    {
        Comparator<Map.Entry<Integer, Integer>> valCmp = new Comparator<Map.Entry<Integer, Integer>>() {
            @Override
            public int compare(Map.Entry<Integer, Integer> o1, Map.Entry<Integer, Integer> o2) {
                return o2.getValue().intValue() - o1.getValue().intValue();  // 升序排序
            }
        };
        List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer, Integer>>(maps.entrySet()); // 传入maps实体
        Collections.sort(list, valCmp);

        return list.get(0).getValue();
    }


    public static int getMapMinValue(Map<Integer, Integer> maps)  // 取出最小value
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

    // 染色体初始化
    private static BitSet generateDataPlacementBitSet(int seversNumber, int minDegree, int degreeDistance) {
        // 生成一个BitSet
        BitSet bitSet = new BitSet(seversNumber);
        // 归一化概率
        mDegrees.forEach((key, value) -> {
            // 每个位置放置数据的归一化概率
            int probability = value - minDegree / degreeDistance;
            // 依照归一化概率放入0或1
            if (Math.random() < probability) bitSet.set(key);

        });
        return bitSet;
    }

    public static BitSet generateChunkNumberBitset(int maxDegree){
        int chunkNumber = new Random().nextInt(maxDegree - 1) + 2;
        String binaryNumber = Integer.toBinaryString(chunkNumber);
        return binaryStringToBitSet(binaryNumber);
    }

    private static BitSet binaryStringToBitSet(String binaryString) {
        BitSet bitSet = new BitSet();

        // 从右到左遍历二进制字符串，将1的位置设置为true
        int length = binaryString.length();
        for (int i = 0; i < length; i++) {
            char bit = binaryString.charAt(length - i - 1);
            if (bit == '1') {
                bitSet.set(i);
            }
        }

        return bitSet;
    }

    // Define the fitness function
    private static double fitness(Genotype<BitGene> PLSolution) {
        Chromosome<BitGene> M_number_chromosome = PLSolution.getChromosome(0);
        Chromosome<BitGene> data_placement_chromosome = PLSolution.getChromosome(1);
        int N_DataNumber = 0;
        int M_chunkNumber = M_number_chromosome.getGene(0).intValue();
        // double smoothingFactor = 1;  // 平滑项，可为其他适当的小数值
        // double penaltyValue = 100; // 惩罚项
        double fitness;
        boolean isFeasibleSolution = true;
        double dataCost;
        int dataBenefit = 0;
        int numberOfSeversLackingData = 0;
        int TotalNumberOfSevers = data_placement_chromosome.length();

        // 校验该数据放置策略是否可行：判断每个服务器是否能访问到足够多的数据块
        // 计算N：总数据块数; 计算数据增益：该数据放置在该位置能服务的服务器数量
        for (int i = 0; i < TotalNumberOfSevers; i++) {
            N_DataNumber += data_placement_chromosome.getGene(i).intValue();
            dataBenefit += data_placement_chromosome.getGene(i).intValue() * mDegrees.get(i);
        }
        if (N_DataNumber < M_chunkNumber) {
            // 如果总数据块不够则直接判定不可行
            isFeasibleSolution = false;
        } else {
            // 总数据块数量足够，开始判断每个服务器需要的数据块是否足够
            for (int i = 0; i < TotalNumberOfSevers; i++) {
                int mRequired = M_chunkNumber;
                for (int j = 0; j < TotalNumberOfSevers; j++) {
                    boolean canAccess = mAccessMatrix[i][j] == 1;
                    boolean hasData = data_placement_chromosome.getGene(j).intValue() == 1;
                    if (canAccess && hasData) --mRequired;
                    if (mRequired == 0) break;
                }
                if (mRequired > 0) {
                    isFeasibleSolution = false;
                    ++numberOfSeversLackingData;
                    break;
                }
                // System.out.println(isFeasibleSolution);
                // System.out.println(mRequired);
            }
        }
        dataCost = -((double) N_DataNumber / M_chunkNumber);
        // 适应度：å成本 + ß是否可行解
        double a = 1;
        double b = 100;
        // double c = 200;
        // System.out.println("isFeasibleSolution  "+ isFeasibleSolution);

        double isFeasibleSolutionValue = (double) numberOfSeversLackingData / TotalNumberOfSevers;
        fitness = a * dataCost + b * isFeasibleSolutionValue;
        if (isFeasibleSolution) {
            System.out.println("Solution： " + M_chunkNumber + ", " + data_placement_chromosome);
        }
        return fitness;
    }

    public void runGACost(int population, int serverNumber) {
        // 首先进行矩阵转化，计算度
        ConvertDistoAccAndCalDegree();
        int maxDegree = getMapMaxValue(mDegrees);
        int minDegree = getMapMinValue(mDegrees);
        int degreeDistance = maxDegree - minDegree;
        BitSet data_place_bitSet = generateDataPlacementBitSet(serverNumber, minDegree, degreeDistance);
        BitSet chunk_number_bitSet = generateChunkNumberBitset(maxDegree);

        BitChromosome data_placement_chromosome = BitChromosome.of(data_place_bitSet);
        BitChromosome M_number_chromosome = BitChromosome.of(chunk_number_bitSet);

        Genotype<BitGene> PLSolution = Genotype.of(M_number_chromosome, data_placement_chromosome);

        System.out.println("PLSolution" + PLSolution);

        Engine<BitGene, Double> engine = Engine.builder(GAModel::fitness, PLSolution)
                .populationSize(population)
                .offspringFraction(0.8)
                .survivorsFraction(0.2)
                .selector(new TournamentSelector<>()) // 选择
                .alterers(new Mutator<>(0.2))  // 变异概率0.2
                .build();

        EvolutionStatistics<Double, ?> Statistics = EvolutionStatistics.ofNumber();

        EvolutionResult<BitGene, Double> result = engine.stream()
                .limit(Limits.bySteadyFitness(50))
                .peek(Statistics)
                .collect(EvolutionResult.toBestEvolutionResult());

        double M = result.getBestPhenotype().getGenotype().getChromosome(0).getGene(0).intValue();
        double N = 0;
        for (int i = 0; i < data_placement_chromosome.length(); i++) {
            N += result.getBestPhenotype().getGenotype().getChromosome(1).getGene(i).intValue();
        }

        System.out.println("Best Solution: " + result.getBestPhenotype());
        System.out.println("GAModel Cost: " + N / M);

    }
}
