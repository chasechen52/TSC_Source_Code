package models;

import java.io.*;
import java.util.*;

//import io.jenetics.EnumChromosome;
import io.jenetics.EnumGene;
import io.jenetics.Genotype;



public class SimpleGAModel {
    private int mServersNumber;   //服务器数量
    private int[][] mAccessMatrix;   //可访问矩阵,初始为全0矩阵
    private int[][] mDistanceMatrix;  //最短距离矩阵
    private int mhops;     //允许的跳数
    private int mPacketsNeed;  //最后需要的分块数

    private double mCost; //Cost
    private double mReplicaCost;
    private Map<Integer,Integer> mDegrees; //度
    private Map<Integer,Integer> mDataPacketsNeed; //尚需数据块的个数
    private List<Integer> mSelectedServerList; //数据节点列表


    public SimpleGAModel(int serversNumber, int [][]distancematrix, int hops) //构造函数
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
    public void initDataPacketsNeed(int requiredpackets)
    {
        for (int key = 0; key < mServersNumber; ++key) //将mDataPacketsNeed重置
        {
            mDataPacketsNeed.put(key,requiredpackets);
        }
    }

    public void ConvertDistoAccAndCalDegree () //根据hops，将能访问到的节点全部置1，否则置0;同时计算度
    {
        for(int i = 0; i < mServersNumber; ++i)
        {
            int access = 0; //可访问节点数
            for(int j = 0; j < mServersNumber; ++j)
            {
                if (mDistanceMatrix[i][j] <= mhops) //当最短距离小于hops时，可访问节点数自增
                {
                    access ++;
                    mAccessMatrix[i][j] = 1;
                }
                else mAccessMatrix[i][j] = 0;
            }
            mDegrees.put(i,access); //access为广义上的度
        }
    }
    public static int getMapMaxValueKey(Map<Integer, Integer> map)  //取出最大value对应的Key值,Int
    {
        List<Map.Entry<Integer,Integer>> list = new ArrayList<Map.Entry<Integer,Integer>>(map.entrySet());
        Collections.sort(list, (o1, o2) -> (o1.getValue().intValue() - o2.getValue().intValue()));
        int key = list.get(list.size() - 1).getKey();
        return key;
    }

    public int getMapMinValue(Map <Integer,Integer> maps)  //取出最小value
    {
        Comparator<Map.Entry<Integer, Integer>> valCmp = new Comparator<Map.Entry<Integer,Integer>>() {
            @Override
            public int compare(Map.Entry<Integer, Integer> o1, Map.Entry<Integer, Integer> o2) {
                return o1.getValue().intValue()-o2.getValue().intValue();  // 升序排序
            }
        };
        List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer,Integer>>(maps.entrySet()); //传入maps实体
        Collections.sort(list,valCmp);

        return list.get(0).getValue();
    }

    public int selectServerwithMostDegrees() //选择有最高度的节点
    {
        return getMapMaxValueKey(mDegrees);
    }

    public boolean checkPacketsRequired() //检查是否所有点都已满足，若是返回true
    {
        int count = 0;
        for (int val:mDataPacketsNeed.values())
        {
            if (val == 0){
                count++;
            }
        }
        if (count == mDataPacketsNeed.size())
        {
            return true;
        }
        return false;
    }

    public static void main(String[] args) {
//        int chromosomeLength = 100; // 染色体长度，即基因的数量
//
//        // 创建枚举染色体，表示包含 100 个基因的染色体，每个基因用 0 或 1 枚举值编码
//        Genotype<EnumGene<BinaryGene>> genotype =
//
//        // 获取染色体
////        EnumChromosome<BinaryGene> enumChromosome = genotype.getChromosome().as(EnumChromosome.class);
//
//        // 输出染色体
//        System.out.println("Enum Chromosome: " + enumChromosome);
//
//        // 或者获取基因的值数组
//        BinaryGene[] geneValues = enumChromosome.toArray(new BinaryGene[0]);
//        System.out.println("Gene Values: " + java.util.Arrays.toString(geneValues));
    }


}
