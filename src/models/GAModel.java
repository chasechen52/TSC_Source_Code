package models;

import io.jenetics.*;
import io.jenetics.engine.*;
import io.jenetics.util.ISeq;
import oracle.jrockit.jfr.events.Bits;

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

    private static int mapMinDegree;

    private static int chunk_chromosome_length;

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
        mapMinDegree = getMapMinValue(mDegrees);
        chunk_chromosome_length = Integer.toBinaryString(mapMinDegree).length();
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

    /**
     * 获取候选服务器List中到已选择服务器List距离最大的TopK个服务器。
     *
     * @param candidateServersList 候选服务器List
     * @param selectedServerList   已选择服务器List
     * @return 候选服务器List中到已选择服务器List距离最大的前TopK个服务器
     */
    public List<Integer> getTopKMaxDistanceServers(List<Integer> candidateServersList, List<Integer> selectedServerList, int TopK) {
        List<Integer> topKMaxDistanceServers = new ArrayList<>(TopK);

        // 创建一个根据最大距离排序的自定义比较器
        Comparator<Integer> maxDistanceComparator = (server1, server2) -> {
            int distance1 = calculateMaxDistance(server1, selectedServerList);
            int distance2 = calculateMaxDistance(server2, selectedServerList);
            return Integer.compare(distance2, distance1); // 降序排列
        };

        // 使用自定义比较器对候选服务器List进行排序
        Collections.sort(candidateServersList, maxDistanceComparator);

        // 从排序后的列表中取前3个服务器
        for (int i = 0; i < Math.min(TopK, candidateServersList.size()); i++) {
            topKMaxDistanceServers.add(candidateServersList.get(i));
        }

        return topKMaxDistanceServers;
    }


    /**
     * 计算候选服务器到已选择服务器List中距离的最大值。
     *
     * @param candidateServer    候选服务器
     * @param selectedServerList 已选择服务器List
     * @return 候选服务器到已选择服务器List中距离的最大值
     */
    private int calculateMaxDistance(int candidateServer, List<Integer> selectedServerList) {
        int maxDistance = Integer.MIN_VALUE;

        for (Integer selectedServer : selectedServerList) {
            int distance = mDistanceMatrix[candidateServer][selectedServer];
            maxDistance = Math.max(maxDistance, distance);
        }

        return maxDistance;
    }

    /**
     * 从topKMaxDistanceServers中随机选择一个数。
     *
     * @param topKMaxDistanceServers topKMaxDistanceServers列表
     * @return 随机选择的一个数
     */
    public int getRandomServerFromTopK(List<Integer> topKMaxDistanceServers) {
        if (topKMaxDistanceServers.isEmpty()) {
            throw new IllegalArgumentException("top3MaxDistanceServers列表不能为空");
        }

        // 使用Random类生成随机数
        Random random = new Random();
        int randomIndex = random.nextInt(topKMaxDistanceServers.size());

        return topKMaxDistanceServers.get(randomIndex);
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
            List<Integer> topKMaxDistanceServers = getTopKMaxDistanceServers(candidateServersList, SelectedServerList, candidatesNumber / 2 + 1); // 获取candidateSeversList中最大的TopK个服务器
            // int initialServer = getLongestDistanceKey(candidateServersList, SelectedServerList); // 在候选服务器列表中选择距离解服务器最大的加入解服务器列表
            int initialServer = getRandomServerFromTopK(topKMaxDistanceServers); // 在候选服务器列表中选择距离解服务器最大的加入解服务器列表
            SelectedServerList.add(initialServer);
            updatePacketsNeed(initialServer); // 选择后更新每个节点所需要的数据包数量
            mDegrees.remove(initialServer); // 如果某节点已经被选择，则在mDegrees中删除
            // updateMDegreesMap(mDegrees, initialServer);
            while (!checkPacketsRequired()) { // 判断当前部署方案是否满足数据请求要求，mDataPacketsNeed是否全为0
                // 更新候选服务器
                candidateServersList = getCandidateServersList(mDegrees, candidatesNumber);
                // System.out.println("candidateServersList" + candidateServersList);
                // System.out.println("candidateServersList : " + candidateServersList);
                List<Integer> newTopKMaxDistanceServers = getTopKMaxDistanceServers(candidateServersList, SelectedServerList, candidatesNumber / 2 + 1); // 获取candidateSeversList中最大的TopK个服务器
                // int newSelectedServer = getLongestDistanceKey(candidateServersList, SelectedServerList);
                int newSelectedServer = getRandomServerFromTopK(newTopKMaxDistanceServers);

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

    public Genotype<BitGene> getSuperiorGenotype() throws IOException, ClassNotFoundException {
        ConvertDistoAccAndCalDegree();

        getReliableSolutions();

        BitSet bitSet = new BitSet(chunk_chromosome_length);
        // System.out.println("chunk_chromosome_length: " + chunk_chromosome_length);
        mSelectedServerList.forEach(index -> bitSet.set(index, true));

        BitChromosome chunk_number_chromosome = BitChromosome.of(Integer.toBinaryString(mSelectedServerList.size()));
        BitChromosome data_placement_chromosome = BitChromosome.of(bitSet);

        return Genotype.of(chunk_number_chromosome, data_placement_chromosome);
    }

    public Genotype<BitGene> getRandomGenotype() throws IOException, ClassNotFoundException {
        ConvertDistoAccAndCalDegree();

        getReliableSolutions();

        BitSet bitSet = new BitSet(chunk_chromosome_length);
        mSelectedServerList.forEach(index -> bitSet.set(index, true));

        BitChromosome chunk_number_chromosome = BitChromosome.of(Integer.toBinaryString(mSelectedServerList.size()));
        BitChromosome data_placement_chromosome = BitChromosome.of(bitSet);

        return Genotype.of(chunk_number_chromosome, data_placement_chromosome);
    }

    public ISeq<Genotype<BitGene>> getInitialPopulation(int populationSize) throws IOException, ClassNotFoundException {
        List<Genotype<BitGene>> genotypeList = new ArrayList<>(populationSize);
        int boundary = populationSize / 2;
        for (int i = 0; i < populationSize; i++) {
            // 一半优秀解，一半随机解
            Genotype<BitGene> genotype = i < boundary ? getSuperiorGenotype() : getRandomGenotype();
            genotypeList.add(genotype);
        }
        return ISeq.of(genotypeList);
    }


    // Define the fitness function
    private static double fitness(Genotype<BitGene> data_placement_strategy) {
        BitChromosome chunk_number_chromosome = data_placement_strategy.getChromosome(0).as(BitChromosome.class);
        BitChromosome data_placement_chromosome = data_placement_strategy.getChromosome(1).as(BitChromosome.class);
        // System.out.println("data_placement_chromosome:   " + data_placement_chromosome);
        // 这里是为了防止比特转换时，出现全为1导致的越界
        // System.out.println("chunk_number_chromosome.intValue(): " + chunk_number_chromosome.intValue());
        // int M = Math.min(chunk_number_chromosome.intValue(), mapMinDegree);
        int M = chromosomeToInt(chunk_number_chromosome);
        System.out.println("M: " + chunk_number_chromosome);
        double N = data_placement_chromosome.bitCount();
        boolean isFeasibleSolution = isFeasibleSolution(chunk_number_chromosome, data_placement_chromosome);
        double fitness;
        // int seversWithoutEnoughData = 0;
        double a = Math.pow(mServersNumber, 2);
        // 适应项：
        double adaptationItem = a * Math.pow(M / N, 2);
        // 惩罚项：不可行的节点数 / 总服务器数
        double penaltyItem = isFeasibleSolution ? 1 : 0;
        // double penaltyItem = b * seversWithoutEnoughData / mServersNumber;
        fitness = adaptationItem * penaltyItem;
        if (isFeasibleSolution) {
            // System.out.println("可行解： " + "M=" +M + ",  数据分布：" + data_placement_chromosome);
        }
        return fitness;
    }

    public static int chromosomeToInt(BitChromosome chromosome) {
        ISeq<BitGene> genes = chromosome.toSeq();
        int value = 0;

        // 将二进制转换为整数
        for (int i = 0; i < genes.length(); i++) {
            if (genes.get(i).getBit()) {
                value |= 1 << (genes.length() - i - 1);
            }
        }

        // 调整范围，确保结果在1到M之间
        // 注意：这里的实现假设解码的值已经是在正确的范围内
        // 如果原始染色体编码方式保证了值的范围，这一步可能不需要或需要根据具体情况调整

        return value + 1; // 如果编码时值域是[1, M]，解码时需要加1
    }

    private static boolean isFeasibleSolution(BitChromosome chunk_number_chromosome, BitChromosome data_placement_chromosome) {
        // 这里是为了防止比特转换时，出现全为1导致的越界
        int M = chromosomeToInt(chunk_number_chromosome);
        // System.out.println("chunk_number_chromosome: " + chunk_number_chromosome);
        // System.out.println("M:  " + M);
        int N = data_placement_chromosome.bitCount();
        boolean isFeasibleSolution = true;
        // 校验该数据放置策略是否可行：判断每个服务器是否能访问到足够多的数据块
        if (N < M) {
            // 如果总数据块不够则直接判定不可行
            isFeasibleSolution = false;
        } else {
            // 总数据块数量足够，开始判断每个服务器需要的数据块是否足够
            for (int i = 0; i < mServersNumber; i++) {
                int mRequired = M;
                for (int j = 0; j < mServersNumber; j++) {
                    boolean canAccess = mAccessMatrix[i][j] == 1;
                    boolean hasData = data_placement_chromosome.getGene(j).getBit();
                    if (canAccess && hasData) --mRequired;
                    if (mRequired == 0) break;
                }
                if (mRequired > 0) {
                    isFeasibleSolution = false;
                    break;
                }
            }
        }
        return isFeasibleSolution;
    }

    public void runGACost(int population, int serverNumber) throws IOException, ClassNotFoundException {
        // 首先进行矩阵转化，计算度
        ConvertDistoAccAndCalDegree();


        // IntegerChromosome data_placement_chromosome = IntegerChromosome.of(0, 1, serverNumber);
        // IntegerChromosome M_number_chromosome = IntegerChromosome.of(2, getMapMinValue(mDegrees), 1);

        BitChromosome chunk_number_chromosome = BitChromosome.of(chunk_chromosome_length);
        BitChromosome data_placement_chromosome = BitChromosome.of(serverNumber);


        Genotype<BitGene> factory = Genotype.of(chunk_number_chromosome, data_placement_chromosome);
        // System.out.println("PLSolution" + PLSolution);

        Engine<BitGene, Double> engine = Engine.builder(GAModel::fitness, factory)
                // .constraint(new SolutionConstraint())
                // .offspringFraction(0.8)
                // .survivorsFraction(0.2)
                // .populationSize(population)
                .populationSize(40)
                // .alterers(new Mutator<>(0.2))  // 变异概率0.2
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

        ISeq<Genotype<BitGene>> initialPopulation = getInitialPopulation(20);

        // System.out.println("initialPopulation: " + initialPopulation.length());

        EvolutionResult<BitGene, Double> result = engine.stream(initialPopulation)
                .limit(Limits.bySteadyFitness(100))
                .peek(Statistics)
                .collect(EvolutionResult.toBestEvolutionResult());

        // System.out.println("population:  " + result.getPopulation().size());

        BitChromosome M_chromosome = result.getBestPhenotype().getGenotype().getChromosome(0).as(BitChromosome.class);
        BitChromosome N_chromosome = result.getBestPhenotype().getGenotype().getChromosome(1).as(BitChromosome.class);

        int M = Math.min(M_chromosome.intValue(), mapMinDegree);
        int N = N_chromosome.bitCount();

        mCost = (double) N / M;

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
