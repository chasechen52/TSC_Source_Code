package models;

import io.jenetics.*;
import io.jenetics.engine.*;
import io.jenetics.util.ISeq;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.concurrent.ConcurrentHashMap;


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

    private static ArrayList<List> mCplexSolutionList;

    private int mPopulation; // 种群数量

    private int mapMinDegree;

    public GAModel(int serversNumber, int[][] distancematrix, int hops) // 构造函数
    {
        mServersNumber = serversNumber;
        mDistanceMatrix = distancematrix;
        mhops = hops;

        mAccessMatrix = new int[mServersNumber][mServersNumber];
        mSelectedServerList = new ArrayList<>();
        mDegrees = new ConcurrentHashMap<>();
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

    private static IntegerGene[] initSolutionGeneArray(List<Integer> solution) {
        IntegerGene[] array = new IntegerGene[mServersNumber];
        IntegerGene geneZero = IntegerGene.of(0, 0, 1);
        IntegerGene geneOne = IntegerGene.of(1, 0, 1);

        Arrays.fill(array, geneZero);
        // 将List中指定索引基因位置设为1
        for (int index : solution) {
            if (index >= 0 && index < mServersNumber) {
                array[index] = geneOne;
            }
        }

        return array;
    }

    private static List<Map.Entry<Integer, Integer>> getRandomHalfEntries(Map<Integer, Integer> map, int topK) {
        List<Map.Entry<Integer, Integer>> shuffledEntries = map.entrySet().stream()
                .collect(Collectors.collectingAndThen(
                        Collectors.toList(),
                        list -> {
                            Collections.shuffle(list);
                            return list.stream();
                        }))
                .collect(Collectors.toList());

        int halfSize = shuffledEntries.size() / 2;

        return shuffledEntries.stream()
                .limit(halfSize < topK ? shuffledEntries.size() : halfSize)
                .collect(Collectors.toList());
    }


    public static List<Integer> getCandidateServersList(Map<Integer, Integer> map, int topK) {
        if (map.isEmpty()) {
            return Collections.emptyList();
        }
        List<Map.Entry<Integer, Integer>> halfEntries = getRandomHalfEntries(map, topK);

        return halfEntries.stream()
                .sorted(Map.Entry.<Integer, Integer>comparingByValue().reversed())
                .limit(Math.min(topK, halfEntries.size()))
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());
    }


    public int getLongestDistanceKey(List<Integer> candidateServersList, List<Integer> selectedServerList) {
        // 计算方法一： 最短距离中的最大距离
        int distanceInSet = Integer.MIN_VALUE;
        int newSelectedServerKey = 0;
        for (Integer candidateServer : candidateServersList) {
            int minDistanceKey = 0;
            int distance = Integer.MAX_VALUE;
            // 计算候选服务器到解服务器的距离
            for (Integer selectedServer : selectedServerList) {
                if (mDistanceMatrix[candidateServer][selectedServer] < distance) {
                    minDistanceKey = selectedServer;
                    distance = mDistanceMatrix[candidateServer][selectedServer];
                }
            }
            if (distanceInSet < distance) {
                distanceInSet = distance;
                newSelectedServerKey = candidateServer;
            }
            // System.out.println("newSelectedServerKey" + newSelectedServerKey);
        }
        return newSelectedServerKey;
    }

    public void updatePacketsNeed(int newSelectedServerKey) {
        // 对newSelectedServerKey可访问节点的mDataPacketsNeed值减1
        for (int server = 0; server < mServersNumber; ++server) {
            if (mAccessMatrix[newSelectedServerKey][server] == 1) {
                int packetsNeed = mDataPacketsNeed.get(server);
                if (packetsNeed > 0) {
                    mDataPacketsNeed.put(server, packetsNeed - 1);
                    // System.out.println("mDataPacketsNeed: " + mDataPacketsNeed);
                }
            }
        }
    }

    public void updateMDegreesMap(Map<Integer, Integer> mDegrees, int selectedServer) {
        // System.out.println("mDegrees" + mDegrees);
        for (int i = 0; i < mServersNumber; ++i) {
            boolean isDeleted = mDegrees.containsKey(selectedServer) && mDegrees.containsKey(i);
            if (isDeleted && mDistanceMatrix[i][selectedServer] <= mhops) // 当最短距离小于hops时，可访问节点数自增
            {
                mDegrees.put(i, mDegrees.get(i) - 1);
            }
            // mDegrees.put(i, access); // access为广义上的度/
        }
    }

    public void getReliableSolutions() throws IOException, ClassNotFoundException {
        // ConvertDistoAccAndCalDegree(); // 首先进行矩阵转化，度计算
        double min_cost = Double.MAX_VALUE;
        int best_pn = 2; // 最低cost对应的数据块数
        int candidatesNumber = 5;
        List<Integer> MinCostServerList = new ArrayList<>();
        // System.out.println("degreesMap" + degreesMap);
        for (int m = 2; m <= mapMinDegree; ++m) {
            mDegrees.clear();
            ConvertDistoAccAndCalDegree(); // 由于每次循环对mDegrees进行了删除，因此每次都需要重新生成mDegrees
            initDataPacketsNeed(m); // 初始化mDataPacketsNeed
            double cost;
            List<Integer> SelectedServerList = new ArrayList<>();
            List<Integer> candidateServersList = getCandidateServersList(mDegrees, candidatesNumber); // step1: 获取候选服务器列表
            int initialServer = getLongestDistanceKey(candidateServersList, SelectedServerList); // 在候选服务器列表中选择距离解服务器最大的加入解服务器列表
            SelectedServerList.add(initialServer);
            updatePacketsNeed(initialServer); // 选择后更新每个节点所需要的数据包数量
            mDegrees.remove(initialServer); // 如果某节点已经被选择，则在mDegrees中删除
            updateMDegreesMap(mDegrees, initialServer);
            while (!checkPacketsRequired()) { // 判断当前部署方案是否满足数据请求要求，mDataPacketsNeed是否全为0
                // 更新候选服务器
                candidateServersList = getCandidateServersList(mDegrees, candidatesNumber);
                // System.out.println("candidateServersList" + candidateServersList);
                // System.out.println("candidateServersList : " + candidateServersList);
                // int newSelectedServer = getLongestDistanceKey(candidateServersList, SelectedServerList);
                int newSelectedServer = getLongestDistanceKey(candidateServersList, SelectedServerList);
                SelectedServerList.add(newSelectedServer);
                updatePacketsNeed(newSelectedServer); // 更新每个服务器需要的数据包数量
                mDegrees.remove(newSelectedServer);  // 如果某节点已经被选择，则在mDegrees中删除
                // updateMDegreesMap(mDegrees, newSelectedServer);
            }
            cost = (double) SelectedServerList.size() / (double) m;
            if (m >= 2 && cost < min_cost) {
                min_cost = cost;
                MinCostServerList = deepCopy(SelectedServerList);
                best_pn = m;
            }
            // System.out.println("SelectedServerList" + SelectedServerList.size());
        }
        // System.out.println("MinCostServerList" + MinCostServerList);
        mSelectedServerList = MinCostServerList;
        mPacketsNeed = best_pn;

        // System.out.println("mSelectedServerList:  " + mSelectedServerList);
        // System.out.println("best_pn: " + mPacketsNeed);
    }

    public Genotype<IntegerGene> getSolutionGenotype() throws IOException, ClassNotFoundException {
        ConvertDistoAccAndCalDegree();
        getReliableSolutions();
        IntegerGene[] geneArray = initSolutionGeneArray(mSelectedServerList);
        // System.out.println("mSelectedServerList: " + mSelectedServerList);
        IntegerChromosome dataPlacement_chromosome = IntegerChromosome.of(geneArray);
        IntegerChromosome M_chromosome = IntegerChromosome.of(2, mapMinDegree, 1);
        return Genotype.of(M_chromosome, dataPlacement_chromosome);
    }

    public ISeq<Genotype<IntegerGene>> getInitialPopulation(int populationSize) throws IOException, ClassNotFoundException {
        Genotype<IntegerGene>[] genotypes = new Genotype[populationSize];
        for (int i = 0; i < populationSize; i++) {
            genotypes[i] = getSolutionGenotype();
            // System.out.println("genotypes: " + genotypes[i].getChromosome(1));
        }
        return ISeq.of(genotypes);
    }


    // Define the fitness function
    private static double fitness(Genotype<IntegerGene> PLSolution) {
        Chromosome<IntegerGene> M_number_chromosome = PLSolution.getChromosome(0);
        Chromosome<IntegerGene> data_placement_chromosome = PLSolution.getChromosome(1);
        // System.out.println("data_placement_chromosome:   " + data_placement_chromosome);
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
        double a = Math.pow(mServersNumber, 4);
        // double b = 10000;
        // 适应项：
        // double adaptationItem = a * N * (1 + (double) 1 / M);
        double adaptationItem = a * Math.pow(M / N, 4);
        // 惩罚项：不可行的节点数 / 总服务器数
        double penaltyItem = isFeasibleSolution ? 1 : 0;
        // double penaltyItem = b * seversWithoutEnoughData / mServersNumber;
        fitness = adaptationItem * penaltyItem;
        if (isFeasibleSolution) {
            // System.out.println("可行解： " + "M=" +M + ",  数据分布：" + data_placement_chromosome);
        }
        return fitness;
    }

    public void runGACost(int population, int serverNumber) throws IOException, ClassNotFoundException {
        // 首先进行矩阵转化，计算度
        ConvertDistoAccAndCalDegree();
        mapMinDegree = getMapMinValue(mDegrees);

        IntegerChromosome data_placement_chromosome = IntegerChromosome.of(0, 1, serverNumber);
        IntegerChromosome M_number_chromosome = IntegerChromosome.of(2, getMapMinValue(mDegrees), 1);
        Genotype<IntegerGene> PLSolution = Genotype.of(M_number_chromosome, data_placement_chromosome);
        // System.out.println("PLSolution" + PLSolution);

        Engine<IntegerGene, Double> engine = Engine.builder(GAModel::fitness, PLSolution)
                // .constraint(new SolutionConstraint())
                // .offspringFraction(0.8)
                // .survivorsFraction(0.2)
                // .populationSize(population)
                .alterers(new Mutator<>(0.2))  // 变异概率0.2
                .offspringFraction(0.3) // 产生子代的比例
                .optimize(Optimize.MAXIMUM) //
                .offspringSelector(new TournamentSelector<>()) // 选择参与交叉变异
                .survivorsSelector(new TournamentSelector<>()) // 选择存活个体
                // .constraint(new Constraint<IntegerGene, Double>() {
                //     @Override
                //     public boolean test(Phenotype<IntegerGene, Double> individual) {
                //         // System.out.println("individual: "  + individual);
                //         Chromosome<IntegerGene> M_number_chromosome = individual.getGenotype().getChromosome(0);
                //         Chromosome<IntegerGene> data_placement_chromosome = individual.getGenotype().getChromosome(1);
                //         double N = 0;
                //         int M = M_number_chromosome.getGene(0).intValue();
                //         boolean isFeasibleSolution = true;
                //
                //         // 校验该数据放置策略是否可行：判断每个服务器是否能访问到足够多的数据块
                //         // 计算N：总数据块数; 计算数据增益：该数据放置在该位置能服务的服务器数量
                //         for (int i = 0; i < mServersNumber; i++) {
                //             N += data_placement_chromosome.getGene(i).intValue();
                //         }
                //         if (N < M) {
                //             // 如果总数据块不够则直接判定不可行
                //             isFeasibleSolution = false;
                //         } else {
                //             // 总数据块数量足够，开始判断每个服务器需要的数据块是否足够
                //             for (int i = 0; i < mServersNumber; i++) {
                //                 int mRequired = M;
                //                 for (int j = 0; j < mServersNumber; j++) {
                //                     boolean canAccess = mAccessMatrix[i][j] == 1;
                //                     boolean hasData = data_placement_chromosome.getGene(j).intValue() == 1;
                //                     if (canAccess && hasData) --mRequired;
                //                     if (mRequired == 0) break;
                //                 }
                //             }
                //         }
                //         return isFeasibleSolution;
                //     }
                //
                //     @Override
                //     public Phenotype<IntegerGene, Double> repair(Phenotype<IntegerGene, Double> individual, long generation) {
                //         boolean isNeedRepair = !test(individual);
                //         if (isNeedRepair) {
                //             // System.out.println(111111111);
                //             try {
                //                 return Phenotype.of(getSolutionGenotype(), generation);
                //             } catch (IOException | ClassNotFoundException e) {
                //                 throw new RuntimeException(e);
                //             }
                //         }
                //         return null;
                //     }
                // })
                .build();

        EvolutionStatistics<Double, ?> Statistics = EvolutionStatistics.ofNumber();

        ISeq<Genotype<IntegerGene>> initialPopulation = getInitialPopulation(30);

        EvolutionResult<IntegerGene, Double> result = engine.stream(initialPopulation)
                .limit(Limits.bySteadyFitness(100))
                .peek(Statistics)
                .collect(EvolutionResult.toBestEvolutionResult());

        double M = result.getBestPhenotype().getGenotype().getChromosome(0).getGene(0).intValue();
        double N = 0;
        for (int i = 0; i < data_placement_chromosome.length(); i++) {
            N += result.getBestPhenotype().getGenotype().getChromosome(1).getGene(i).intValue();
        }

        mCost = N / M;

        // System.out.println("statistics \n" + Statistics);
        // System.out.println("Best Solution: " + result.getBestPhenotype());
        // System.out.println("Best Solution Fitness: " + result.getBestPhenotype().getFitness());
        // System.out.println("GAModel Cost: " + N/M);

    }


    public double getCost() {
        return mCost;
    }

    public int getPacketsNeed() {
        return mPacketsNeed;
    }
}
