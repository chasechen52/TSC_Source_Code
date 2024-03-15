package models;

import org.jgap.*;
import org.jgap.impl.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;


public class JGAPGAModel {

    private static int mServersNumber;   // 服务器数量
    private static int[][] mAccessMatrix;   // 可访问矩阵,初始为全0矩阵
    private int[][] mDistanceMatrix;  // 最短距离矩阵
    private int mhops;     // 允许的跳数
    private int mPacketsNeed;  // 最后需要的分块数

    private double mCost; // Cost
    private double mReplicaCost;
    private static ConcurrentHashMap<Integer, Integer> mDegrees; // 度
    private Map<Integer, Integer> mDataPacketsNeed; // 尚需数据块的个数
    private List<Integer> mSelectedServerList; // 数据节点列表

    private static ArrayList<List> mCplexSolutionList;

    private int mPopulation; // 种群数量

    private int mapMinDegree;

    private double mFitness;

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

        // 报错的原因为此处list为【】
        // System.out.println("list" + list);
        return list.get(0).getValue();
    }

    public JGAPGAModel(int serversNumber, int[][] distancematrix, int hops) // 构造函数
    {
        mServersNumber = serversNumber;
        mDistanceMatrix = distancematrix;
        mhops = hops;

        mAccessMatrix = new int[mServersNumber][mServersNumber];
        mSelectedServerList = new ArrayList<>();
        mDegrees = new ConcurrentHashMap<>();
        mDataPacketsNeed = new HashMap<>();
    }

    public int selectServerwithMostDegrees() // 选择有最高度的节点
    {
        return getMapMaxValueKey(mDegrees);
    }

    // 随机选择一半元素
    private static List<Map.Entry<Integer, Integer>> getRandomHalfEntries(Map<Integer, Integer> map) {
        List<Map.Entry<Integer, Integer>> shuffledEntries = map.entrySet().stream()
                .collect(Collectors.collectingAndThen(
                        Collectors.toList(),
                        list -> {
                            Collections.shuffle(list);
                            return list.stream();
                        }))
                .collect(Collectors.toList());

        int halfSize = shuffledEntries.size() / 2 + 1;

        // System.out.println("halfSize " + halfSize);

        return shuffledEntries.stream().limit(halfSize).collect(Collectors.toList());
    }


    // 获取候选服务器列表，选择最大的topK个
    private static List<Integer> getCandidateServersList(Map<Integer, Integer> map, int topK) {
        if (map.isEmpty()) {
            return Collections.emptyList();
        }
        List<Map.Entry<Integer, Integer>> halfEntries = getRandomHalfEntries(map);

        return halfEntries.stream()
                .sorted(Map.Entry.<Integer, Integer>comparingByValue().reversed())
                .limit(Math.min(topK, halfEntries.size()))
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());
    }


    // 从List中随机选择一个值
    private static int getRandomKeyFromCandidateServersList(List<Integer> valuesList) {
        // System.out.println("valuesList: " +valuesList );
        // 判断List是否为空
        if (valuesList == null || valuesList.isEmpty()) {
            throw new IllegalArgumentException("List不能为空");
        }

        // 生成一个随机索引
        Random random = new Random();
        int randomIndex = random.nextInt(valuesList.size());

        // 返回对应的随机值
        return valuesList.get(randomIndex);
    }

    public static int getMapMaxValueKey(Map<Integer, Integer> map)  // 取出最大value对应的Key值,Int
    {
        List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer, Integer>>(map.entrySet());
        Collections.sort(list, (o1, o2) -> (o1.getValue().intValue() - o2.getValue().intValue()));
        int key = list.get(list.size() - 1).getKey();
        return key;
    }

    public boolean checkPacketsRequired() // 检查是否所有点都已满足，若是返回true
    {
        int count = 0;
        for (int val : mDataPacketsNeed.values()) {
            // System.out.println("mDataPacketsNeed" + mDataPacketsNeed);
            if (val == 0) {
                count++;
            }
        }
        if (count == mDataPacketsNeed.size()) {
            return true;
        }
        return false;
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

            double a = Math.pow(mServersNumber, 2);
            // double b = 10000;
            // 适应项：
            // double adaptationItem = a * N * (1 + (double) 1 / M);
            double adaptationItem = a * Math.pow(M / N, 2);
            // 惩罚项：不可行的节点数 / 总服务器数
            double penaltyItem = isFeasibleSolution ? 1 : 0;
            // double penaltyItem = b * seversWithoutEnoughData / mServersNumber;
            fitness = adaptationItem * penaltyItem;
            if (isFeasibleSolution) {
                ArrayList<Integer> GAServerList = new ArrayList<>();
                for (int i = 0; i < strategy_chromosome.getGenes().length; i++) {
                    Gene gene = strategy_chromosome.getGene(i);
                    int GAServer = (int) gene.getAllele();
                    GAServerList.add(GAServer);
                }
                // System.out.println("can solve GAServerList" + GAServerList);
            }
            return fitness;
            // return 0;
        }
    }


    public void initDataPacketsNeed(int requiredPackets) {
        for (int key = 0; key < mServersNumber; ++key) // 将mDataPacketsNeed重置
        {
            mDataPacketsNeed.put(key, requiredPackets);
        }
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
        for (int server = 0; server < mServersNumber; ++server) {  // 对newSelectedServerKey可访问节点的mDataPacketsNeed值减1
            if (mAccessMatrix[newSelectedServerKey][server] == 1) {
                int packetsNeed = mDataPacketsNeed.get(server);
                if (packetsNeed > 0) {
                    mDataPacketsNeed.put(server, packetsNeed - 1);
                }
            }
        }
    }

    public void updateMDegreesMap(Map<Integer, Integer> mDegrees, int selectedServer) {
        for (int i = 0; i < mServersNumber; ++i) {
            boolean isDeleted = mDegrees.containsKey(selectedServer);
            if (isDeleted && mDistanceMatrix[i][selectedServer] <= mhops) // 当最短距离小于hops时，可访问节点数自增
            {
                mDegrees.put(i, mDegrees.get(i) - 1);
            }
            // mDegrees.put(i, access); // access为广义上的度/
        }
    }

    public void getReliableSolutions() throws IOException, ClassNotFoundException {
        ConvertDistoAccAndCalDegree(); // 首先进行矩阵转化，度计算
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
            // updateMDegreesMap(mDegrees, initialServer);
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
                // updateMDegreesMap(mDegrees,newSelectedServer);
            }
            cost = (double) SelectedServerList.size() / (double) m;
            if (m >= 2 && cost < min_cost) {
                min_cost = cost;
                MinCostServerList = deepCopy(SelectedServerList);
                best_pn = m;
            }
            // System.out.println("SelectedServerList" + SelectedServerList);
        }
        System.out.println("MinCostServerList" + MinCostServerList.size());
        mSelectedServerList = MinCostServerList;
        mPacketsNeed = best_pn;
        // System.out.println("best_pn: " + mPacketsNeed);
    }

    public IChromosome convertSolutionToChromosome(Configuration configuration) throws InvalidConfigurationException {

        Gene[] genes = new Gene[mServersNumber + 1];
        // TODO: mDegrees 删除的问题
        genes[0] = new IntegerGene(configuration, 2, mapMinDegree);
        genes[0].setAllele(mPacketsNeed);
        // System.out.println("genes[j]" + genes[0]);
        for (int j = 1; j < mServersNumber + 1; j++) {
            genes[j] = new IntegerGene(configuration, 0, 1);
            genes[j].setAllele(0);
        }

        for (int mIndex : mSelectedServerList) {
            // int mIndex = (int) index;
            if (mIndex >= 0 && mIndex < mServersNumber + 1) {
                // genes[mIndex] = new IntegerGene(configuration, 0, 1);
                genes[mIndex].setAllele(1);
            }
        }
        return new Chromosome(configuration, genes);
    }


    public void runGACost(int populationSize, int serverNumber) throws InvalidConfigurationException, IOException, ClassNotFoundException {
        Configuration configuration = new DefaultConfiguration();

        configuration.setPopulationSize(populationSize);

        // 设置适应度函数
        FitnessFunction fitnessFunction = new MyFitnessFunction();
        configuration.setFitnessFunction(fitnessFunction);

        // 设置演化方式为精英选择
        WeightedRouletteSelector selectionMethod = new WeightedRouletteSelector(configuration);
        configuration.addNaturalSelector(selectionMethod, false);

        // 设置变异概率为0.2
        int mutationRate = 10;
        MutationOperator mutationOperator = new MutationOperator(configuration, mutationRate);
        configuration.addGeneticOperator(mutationOperator);

        Gene[] genes = new Gene[serverNumber + 1];

        ConvertDistoAccAndCalDegree();
        mapMinDegree = getMapMinValue(mDegrees);

        // 创建复合基因并将两个基因合并成一个染色体
        for (int i = 0; i < genes.length; i++) {
            if (i == 0) genes[i] = new IntegerGene(configuration, 2, mapMinDegree);
            else genes[i] = new IntegerGene(configuration, 0, 1);
        }

        // 创建染色体并将多个基因添加到其中
        IChromosome data_placement_strategy_chromosome = new Chromosome(configuration, genes);

        configuration.setSampleChromosome(data_placement_strategy_chromosome);
        configuration.setPopulationSize(populationSize);

        // 使用Genotype.randomInitialGenotype()创建一个随机初始化的个体的基因型
        Genotype genotype = Genotype.randomInitialGenotype(configuration);
        Population population = genotype.getPopulation();
        for (int i = 0; i < 15; i++) {
            getReliableSolutions();
            ConvertDistoAccAndCalDegree(); // 由于每次循环对mDegrees进行了删除，因此每次都需要重新生成mDegrees
            IChromosome feasible_solution = convertSolutionToChromosome(configuration);
            population.addChromosome(feasible_solution);
        }

        // 执行一定数量的进化操作
        int numberOfEvolutions = 200; // 你可以根据需要设置进化的次数
        for (int i = 0; i < numberOfEvolutions; i++) {
            genotype.evolve();
        }
        // 获取最优染色体
        IChromosome fittestChromosome = genotype.getFittestChromosome();
        double fitness = genotype.getFittestChromosome().getFitnessValue();
        mFitness = fitness;

        double cost;
        double N = 0;

        for (int i = 1; i < serverNumber + 1; i++) {
            N += (int) fittestChromosome.getGene(i).getAllele();
        }
        int M = (int) fittestChromosome.getGene(0).getAllele();

        cost = N / M;

        mCost = cost;

        // 打印最优染色体的基因值
        // System.out.println("mGACost: " + cost);
        // System.out.println("fitness: " + fitness);
        // System.out.println("mGAServers: " + Arrays.toString(fittestChromosome.));
        ArrayList<Integer> GAServerList = new ArrayList<>();
        for (int i = 0; i < fittestChromosome.getGenes().length; i++) {
            Gene gene = fittestChromosome.getGene(i);
            int GAServer = (int) gene.getAllele();
            GAServerList.add(GAServer);
        }

        Configuration.reset();
        // System.out.println("m:  " + fittestChromosome.getGene(0));
        // System.out.println("GAServerList" + GAServerList);
    }


    public double getCost() {
        return mCost;
    }

    public double getFitness() {
        return mFitness;
    }
}







