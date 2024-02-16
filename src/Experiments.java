import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Map.Entry;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import java.text.NumberFormat;

import models.*;
import org.jgap.InvalidConfigurationException;
import tool.RandomGraphGenerator;

public class Experiments {

    private static ECGreedyVoteModel mECGreedyVoteModel;
    private static ECGreedyDegreeModel mECGreedyDegreeModel;
    private static ECRandomModel mECRandomModel;
    private static ECRandomModel mReplicaModel;

    private static GAModel mGAModel;

    private static JGAPGAModel mJGAPGAModel;
    private static double mGACost;
    private static List<Integer> mGAServers = new ArrayList<>();
    private static double mECGreedyVoteCost;
    private static List<Integer> mECGreedyVoteServers = new ArrayList<>();

    private static List<ECCplexModel> mECCplexModels = new ArrayList<>();
    private static double mECCplexCost;
    private static List<Integer> mECCplexServers = new ArrayList<>();

    private static double mECGreedyDegreeCost;
    private static List<Integer> mECGreedyDegreeServers = new ArrayList<>();

    private static double mECRandomCost;
    private static List<Integer> mECRandomServers = new ArrayList<>();

    private static double mReplicaCost;
    private static List<Integer> mReplicaServers = new ArrayList<>();


    private static List<String> mLines = new ArrayList<>();

    private static ArrayList<List> CplexSolutionList = new ArrayList<>();
    private static ArrayList<Integer> voteSolutionList = new ArrayList<>();

    private static double GAFitness;


    public static void main(String[] args) throws ClassNotFoundException, IOException, InvalidConfigurationException {

        // runSettings();
        // runExample();
        // runGAExample();
        runJGAPGAExample();
        try {
            Calendar cal = Calendar.getInstance();
            SimpleDateFormat sdf = new SimpleDateFormat("dd-HH-mm");

            Path file = Paths.get("ECModels Comparision" + sdf.format(cal.getTime()) + ".txt");
            // Path file = Paths.get("Synthetic - results-" + sdf.format(cal.getTime()) + ".txt");

            Files.write(file, mLines, Charset.forName("UTF-8"));
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private static boolean isConnected(List<Integer> servers, int server, int[][] adjacencyMatrix) {
        for (int s : servers) {
            if (adjacencyMatrix[s][server] == 1)
                return true;
        }

        return false;
    }

    private static int getMapMinValue(Map<Integer, Integer> maps)  // 取出最小value对应的Key值,Int
    {
        Comparator<Entry<Integer, Integer>> valCmp = new Comparator<Entry<Integer, Integer>>() {
            @Override
            public int compare(Entry<Integer, Integer> o1, Entry<Integer, Integer> o2) {
                return o1.getValue().intValue() - o2.getValue().intValue();  // 升序排序
            }
        };
        List<Entry<Integer, Integer>> list = new ArrayList<Entry<Integer, Integer>>(maps.entrySet()); // 传入maps实体
        Collections.sort(list, valCmp);

        return list.get(0).getValue();
    }

    private static int minDegree(int ServersNumber, int[][] distanceMatrix, int hops) // 根据hops，将能访问到的节点全部置1，否则置0;同时计算度
    {
        Map<Integer, Integer> Degrees = new HashMap<>();
        for (int i = 0; i < ServersNumber; ++i) {
            int access = 0; // 可访问节点数
            for (int j = 0; j < ServersNumber; ++j) {
                if (distanceMatrix[i][j] <= hops) // 当最短距离小于hops时，可访问节点数自增
                {
                    access++;
                }
            }
            Degrees.put(i, access); // access为广义上的度
        }
        int mindeg = getMapMinValue(Degrees);
        return mindeg;
    }

    private static int[][] GraphGenerate(int serversNumber, double density) {
        RandomGraphGenerator graphGenerator = new RandomGraphGenerator(serversNumber, density);
        graphGenerator.createRandomGraph();
        int[][] adjacencyMatrix = graphGenerator.getRandomGraphAdjacencyMatrix();
        int[][] distanceMatrics = graphGenerator.getRandomGraphDistanceMatrix();

        // int[][] adjacencyMatrix = {
        //         {0, 1, 0, 0, 0, 1, 0, 0, 1, 0},
        //         {1, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        //         {0, 0, 0, 1, 0, 0, 1, 1, 0, 0},
        //         {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        //         {0, 0, 0, 0, 0, 0, 0, 1, 1, 0},
        //         {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        //         {0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        //         {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        //         {1, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        //         {0, 1, 0, 0, 0, 0, 1, 0, 0, 0},
        // };
        //
        // int[][] distanceMatrics = {
        //         {0, 1, 4, 5, 2, 1, 3, 3, 1, 2},
        //         {1, 0, 3, 4, 3, 2, 2, 4, 2, 1},
        //         {4, 3, 0, 1, 2, 5, 1, 1, 3, 2},
        //         {5, 4, 1, 0, 3, 6, 2, 2, 4, 3},
        //         {2, 3, 2, 3, 0, 3, 3, 1, 1, 4},
        //         {1, 2, 5, 6, 3, 0, 4, 4, 2, 3},
        //         {3, 2, 1, 2, 3, 4, 0, 2, 4, 1},
        //         {3, 4, 1, 2, 1, 4, 2, 0, 2, 3},
        //         {1, 2, 3, 4, 1, 2, 4, 2, 0, 3},
        //         {2, 1, 2, 3, 4, 3, 1, 3, 3, 0},
        // };


        // System.out.println("\n------- AdjMatrix -------"); // 打印邻接矩阵方便绘图
        // for (int i = 0; i < serversNumber; ++i) {
        //     System.out.println(Arrays.toString(adjacencyMatrix[i]) + ",");
        // }
        // System.out.println("\n");
        //
        // System.out.println("\n------- distanceMatrics -------"); // 打印邻接矩阵方便绘图
        // for (int i = 0; i < serversNumber; ++i) {
        //     System.out.println(Arrays.toString(distanceMatrics[i]) + ",");
        // }
        // System.out.println("\n");

        return distanceMatrics;
    }

    private static void ModelSetup(int serversNumber, int[][] distanceMatrics, int hops) {
        mECGreedyVoteModel = new ECGreedyVoteModel(serversNumber, distanceMatrics, hops);
        mECGreedyDegreeModel = new ECGreedyDegreeModel(serversNumber, distanceMatrics, hops);
        mECRandomModel = new ECRandomModel(serversNumber, distanceMatrics, hops);
        mReplicaModel = new ECRandomModel(serversNumber, distanceMatrics, hops);
        mGAModel = new GAModel(serversNumber, distanceMatrics, hops);
        mJGAPGAModel = new JGAPGAModel(serversNumber, distanceMatrics, hops);
        int mindeg = minDegree(serversNumber, distanceMatrics, hops);
        // System.out.println("mindeg" + mindeg);
        for (int m = 1; m <= mindeg; m++) {
            ECCplexModel tmp = new ECCplexModel(serversNumber, distanceMatrics, hops, m);
            mECCplexModels.add(tmp);
        }
    }

    private static void runExample() throws ClassNotFoundException, IOException {
        int s = 10;
        double d = 2;
        int h = 3;
        while (true) {
            int[][] dism = GraphGenerate(s, d);
            ModelSetup(s, dism, h);
            runECCplexCost();
            if (mECCplexCost <= mReplicaCost * 0.75) {
                System.out.println("Replica:" + mReplicaCost + "  " + mReplicaServers);
                System.out.println("Erasure Code:" + mECCplexCost + "  " + mECCplexServers);
                break;
            }
        }
    }

    private static void runSettings() throws ClassNotFoundException, IOException {
        // in this set, random userNumber, same density, different numbers of edge
        // server
        int avgepoch = 10;
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(3);  // 小数点后最多保留3位

        int s0 = 150;
        double d0 = 2.0;
        int h0 = 2;
        System.out.println("\nServernumber time test");
        for (int s = 50; s <= 250; s += 50) {
            System.out.println("\nServernumber = " + s + " -------------");
            int[][] dism = GraphGenerate(s, d0);
            ModelSetup(s, dism, h0);

            double replicatime = 0;
            // double cplextime = 0;
            double votetime = 0;
            double degreetime = 0;
            double randomtime = 0;

            long start1 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runReplicaCost();
            }
            long end1 = System.currentTimeMillis();
            replicatime = (end1 - start1) / (double) 10;


            // long start2 = System.currentTimeMillis();
            // for(int i = 0; i < avgepoch; ++i)
            // {
            // 	runECCplexCost();
            // }
            // long end2 = System.currentTimeMillis();
            // cplextime = (end2-start2)/(double) 10;

            long start3 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runECGreedyVoteCost();
            }
            long end3 = System.currentTimeMillis();
            votetime = (end3 - start3) / (double) 10;

            long start4 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runECGreedyDegreeCost();
            }
            long end4 = System.currentTimeMillis();
            degreetime = (end4 - start4) / (double) 10;

            long start5 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runECRandomCost();
            }
            long end5 = System.currentTimeMillis();
            randomtime = (end5 - start5) / (double) 10;

            System.out.println("Time: " + replicatime + " | " + votetime + ", " + degreetime + ", " + randomtime);
        }
        System.out.println("\n");


        System.out.println("Dentisy time test");
        for (double d = 2.0; d <= 5.0; d += 0.6) {
            System.out.println("\nDentisy = " + d + " -------------");
            int[][] dism = GraphGenerate(s0, d);
            ModelSetup(s0, dism, h0);
            double replicatime = 0;
            // double cplextime = 0;
            double votetime = 0;
            double degreetime = 0;
            double randomtime = 0;

            long start1 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runReplicaCost();
            }
            long end1 = System.currentTimeMillis();
            replicatime = (end1 - start1) / (double) 10;


            // long start2 = System.currentTimeMillis();
            // for(int i = 0; i < avgepoch; ++i)
            // {
            // 	runECCplexCost();
            // }
            // long end2 = System.currentTimeMillis();
            // cplextime = (end2-start2)/(double) 10;

            long start3 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runECGreedyVoteCost();
            }
            long end3 = System.currentTimeMillis();
            votetime = (end3 - start3) / (double) 10;

            long start4 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runECGreedyDegreeCost();
            }
            long end4 = System.currentTimeMillis();
            degreetime = (end4 - start4) / (double) 10;

            long start5 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runECRandomCost();
            }
            long end5 = System.currentTimeMillis();
            randomtime = (end5 - start5) / (double) 10;

            System.out.println("Time: " + replicatime + " | " + votetime + ", " + degreetime + ", " + randomtime);
        }
        System.out.println("\n");


        System.out.println("Hops time test");
        for (int h = 1; h <= 5; h++) {
            System.out.println("\nHops = " + h + " -------------");
            int[][] dism = GraphGenerate(s0, d0);
            ModelSetup(s0, dism, h);

            double replicatime = 0;
            // double cplextime = 0;
            double votetime = 0;
            double degreetime = 0;
            double randomtime = 0;

            long start1 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runReplicaCost();
            }
            long end1 = System.currentTimeMillis();
            replicatime = (end1 - start1) / (double) 10;


            // long start2 = System.currentTimeMillis();
            // for(int i = 0; i < avgepoch; ++i)
            // {
            // 	runECCplexCost();
            // }
            // long end2 = System.currentTimeMillis();
            // cplextime = (end2-start2)/(double) 10;

            long start3 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runECGreedyVoteCost();
            }
            long end3 = System.currentTimeMillis();
            votetime = (end3 - start3) / (double) 10;

            long start4 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runECGreedyDegreeCost();
            }
            long end4 = System.currentTimeMillis();
            degreetime = (end4 - start4) / (double) 10;

            long start5 = System.currentTimeMillis();
            for (int i = 0; i < avgepoch; ++i) {
                runECRandomCost();
            }
            long end5 = System.currentTimeMillis();
            randomtime = (end5 - start5) / (double) 10;

            System.out.println("Time: " + replicatime + " | " + votetime + ", " + degreetime + ", " + randomtime);
        }


        System.out.println("-------------- Setting 1 --------------");
        mLines.add("-------------- Setting 1 --------------");
        double d1 = 1.0;
        int s1 = 20;
        int h1 = 2;
        int[][] hopdism1 = GraphGenerate(s1, d1); // 测试hops时保证同一张图
        System.out.println("----Servernumber Test 1----");
        mLines.add("----Servernumber Test 1----");
        for (int s = 10; s <= 35; s += 5) {
            double replicacosts = 0.0;
            double cplexcosts = 0.0;
            double votecosts = 0.0;
            double degreecosts = 0.0;
            double randomcosts = 0.0;
            for (int i = 0; i < avgepoch; ++i) {
                int[][] dism = GraphGenerate(s, d1);
                ModelSetup(s, dism, h1);

                runECCplexCost();
                runECGreedyVoteCost();
                runECGreedyDegreeCost();

                cplexcosts += mECCplexCost;
                votecosts += mECGreedyVoteCost;
                degreecosts += mECGreedyDegreeCost;


                double randcurcost = 0.0;
                double repcurcost = 0.0;
                for (int j = 0; j < avgepoch; ++j) {
                    runReplicaCost();
                    runECRandomCost();
                    randcurcost += mECRandomCost;
                    repcurcost += mReplicaCost;
                }
                randcurcost /= (double) avgepoch;
                repcurcost /= (double) avgepoch;

                replicacosts += repcurcost;
                randomcosts += randcurcost;
            }

            replicacosts /= (double) avgepoch;
            cplexcosts /= (double) avgepoch;
            votecosts /= (double) avgepoch;
            degreecosts /= (double) avgepoch;
            randomcosts /= (double) avgepoch;

            System.out.println("ServersNumber:" + s + "  Density:" + nf.format(d1) + "  Hops:" + h1 + "  [replicacosts | cplexcosts, votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(cplexcosts) + ", " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
            mLines.add("ServersNumber:" + s + "  Density:" + nf.format(d1) + "  Hops:" + h1 + "  [replicacosts | cplexcosts, votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(cplexcosts) + ", " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
        }
        mLines.add("\n");
        System.out.println("\n");

        System.out.println("----Denstiy Test 1----");
        mLines.add("----Denstiy Test 1----");
        for (double d = 1.0; d <= 2.5; d += 0.3) {
            double replicacosts = 0.0;
            double cplexcosts = 0.0;
            double votecosts = 0.0;
            double degreecosts = 0.0;
            double randomcosts = 0.0;
            for (int i = 0; i < avgepoch; ++i) {
                int[][] dism = GraphGenerate(s1, d);
                ModelSetup(s1, dism, h1);

                runECCplexCost();
                runECGreedyVoteCost();
                runECGreedyDegreeCost();

                cplexcosts += mECCplexCost;
                votecosts += mECGreedyVoteCost;
                degreecosts += mECGreedyDegreeCost;


                double randcurcost = 0.0;
                double repcurcost = 0.0;
                for (int j = 0; j < avgepoch; ++j) {
                    runReplicaCost();
                    runECRandomCost();
                    randcurcost += mECRandomCost;
                    repcurcost += mReplicaCost;
                }
                randcurcost /= (double) avgepoch;
                repcurcost /= (double) avgepoch;

                replicacosts += repcurcost;
                randomcosts += randcurcost;

            }

            replicacosts /= (double) avgepoch;
            cplexcosts /= (double) avgepoch;
            votecosts /= (double) avgepoch;
            degreecosts /= (double) avgepoch;
            randomcosts /= (double) avgepoch;

            System.out.println("ServersNumber:" + s1 + "  Density:" + nf.format(d) + "  Hops:" + h1 + "  [replicacosts | cplexcosts, votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(cplexcosts) + ", " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
            mLines.add("ServersNumber:" + s1 + "  Density:" + nf.format(d) + "  Hops:" + h1 + "  [replicacosts | cplexcosts, votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(cplexcosts) + ", " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
        }
        mLines.add("\n");
        System.out.println("\n");

        System.out.println("----Hops Test 1----");
        mLines.add("----Hops Test 1----");

        for (int h = 1; h <= 5; ++h) {
            double replicacosts = 0.0;
            double cplexcosts = 0.0;
            double votecosts = 0.0;
            double degreecosts = 0.0;
            double randomcosts = 0.0;
            ModelSetup(s1, hopdism1, h);


            runECCplexCost();
            runECGreedyVoteCost();
            runECGreedyDegreeCost();

            cplexcosts = mECCplexCost;
            votecosts = mECGreedyVoteCost;
            degreecosts = mECGreedyDegreeCost;

            for (int i = 0; i < avgepoch; ++i) {
                runECRandomCost();
                runReplicaCost();
                replicacosts += mReplicaCost;
                randomcosts += mECRandomCost;
            }
            randomcosts /= (double) avgepoch;
            replicacosts /= (double) avgepoch;

            System.out.println("ServersNumber:" + s1 + "  Density:" + nf.format(d1) + "  Hops:" + h + "  [replicacosts | cplexcosts, votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(cplexcosts) + ", " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
            mLines.add("ServersNumber:" + s1 + "  Density:" + nf.format(d1) + "  Hops:" + h + "  [replicacosts | cplexcosts, votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(cplexcosts) + ", " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
        }
        mLines.add("\n");
        System.out.println("\n");


        mLines.add("\n");
        System.out.println("\n");
        System.out.println("-------------- Setting 2 --------------");
        mLines.add("-------------- Setting 2 --------------");
        double d2 = 2.0;
        int s2 = 150;
        int h2 = 2;
        int[][] hopdism2 = GraphGenerate(s2, d2); // 测试hops时保证同一张图
        System.out.println("----Servernumber Test 2----");
        mLines.add("----Servernumber Test 2----");
        for (int s = 50; s <= 250; s += 50) {
            double replicacosts = 0.0;
            double votecosts = 0.0;
            double degreecosts = 0.0;
            double randomcosts = 0.0;
            for (int i = 0; i < avgepoch; ++i) {
                int[][] dism = GraphGenerate(s, d2);
                ModelSetup(s, dism, h2);

                runECGreedyVoteCost();
                runECGreedyDegreeCost();
                votecosts += mECGreedyVoteCost;
                degreecosts += mECGreedyDegreeCost;


                double randcurcost = 0.0;
                double repcurcost = 0.0;
                for (int j = 0; j < avgepoch; ++j) {
                    runReplicaCost();
                    runECRandomCost();
                    randcurcost += mECRandomCost;
                    repcurcost += mReplicaCost;
                }
                randcurcost /= (double) avgepoch;
                repcurcost /= (double) avgepoch;

                replicacosts += repcurcost;
                randomcosts += randcurcost;

            }

            replicacosts /= (double) avgepoch;
            votecosts /= (double) avgepoch;
            degreecosts /= (double) avgepoch;
            randomcosts /= (double) avgepoch;

            System.out.println("ServersNumber:" + s + "  Density:" + nf.format(d2) + "  Hops:" + h2 + "  [replicacosts | votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
            mLines.add("ServersNumber:" + s + "  Density:" + nf.format(d2) + "  Hops:" + h2 + "  [replicacosts | votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
        }
        mLines.add("\n");
        System.out.println("\n");

        System.out.println("----Denstiy Test 2----");
        mLines.add("----Denstiy Test 2----");
        for (double d = 2.0; d <= 5.0; d += 0.6) {
            double replicacosts = 0.0;
            double votecosts = 0.0;
            double degreecosts = 0.0;
            double randomcosts = 0.0;
            for (int i = 0; i < avgepoch; ++i) {
                int[][] dism = GraphGenerate(s2, d);
                ModelSetup(s2, dism, h2);

                runECGreedyVoteCost();
                runECGreedyDegreeCost();

                votecosts += mECGreedyVoteCost;
                degreecosts += mECGreedyDegreeCost;


                double randcurcost = 0.0;
                double repcurcost = 0.0;
                for (int j = 0; j < avgepoch; ++j) {
                    runReplicaCost();
                    runECRandomCost();
                    randcurcost += mECRandomCost;
                    repcurcost += mReplicaCost;
                }
                randcurcost /= (double) avgepoch;
                repcurcost /= (double) avgepoch;

                replicacosts += repcurcost;
                randomcosts += randcurcost;

            }

            replicacosts /= (double) avgepoch;
            votecosts /= (double) avgepoch;
            degreecosts /= (double) avgepoch;
            randomcosts /= (double) avgepoch;

            System.out.println("ServersNumber:" + s2 + "  Density:" + nf.format(d) + "  Hops:" + h2 + "  [replicacosts | votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
            mLines.add("ServersNumber:" + s2 + "  Density:" + nf.format(d) + "  Hops:" + h2 + "  [replicacosts | votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
        }
        mLines.add("\n");
        System.out.println("\n");

        System.out.println("----Hops Test 2----");
        mLines.add("----Hops Test 2----");

        for (int h = 1; h <= 5; ++h) {
            double replicacosts = 0.0;
            double votecosts = 0.0;
            double degreecosts = 0.0;
            double randomcosts = 0.0;

            ModelSetup(s2, hopdism2, h);

            runECGreedyVoteCost();
            runECGreedyDegreeCost();

            votecosts = mECGreedyVoteCost;
            degreecosts = mECGreedyDegreeCost;

            for (int i = 0; i < avgepoch; ++i) {
                runECRandomCost();
                runReplicaCost();
                replicacosts += mReplicaCost;
                randomcosts += mECRandomCost;
            }
            randomcosts /= (double) avgepoch;
            replicacosts /= (double) avgepoch;


            System.out.println("ServersNumber:" + s2 + "  Density:" + nf.format(d2) + "  Hops:" + h + "  [replicacosts | votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
            mLines.add("ServersNumber:" + s2 + "  Density:" + nf.format(d2) + "  Hops:" + h + "  [replicacosts | votecosts, degreecosts, randomcosts] : [" + nf.format(replicacosts) + " | " + nf.format(votecosts) + ", " + nf.format(degreecosts) + ", " + nf.format(randomcosts) + "]");
        }

    }

    // private static void runGAExample() throws ClassNotFoundException, IOException {
    //     int serverNumber = 20;
    //     double d = 2;
    //     int h = 1;
    //     int population = 100;
    //     while (true) {
    //         int[][] dism = GraphGenerate(serverNumber, d);
    //         ModelSetup(serverNumber, dism, h);
    //         getECCplexSolution();
    //         System.out.println("CplexSolutionList    " + CplexSolutionList);
    //         runGACost(population, serverNumber, CplexSolutionList);
    //         // runJGAPGACost(population, serverNumber, CplexSolutionList);
    //         runECGreedyVoteCost();
    //         // runECCplexCost();
    //         break;
    //         // if (mGACost <= mReplicaCost * 0.75) {
    //         //     System.out.println("Erasure Code Cplex:" + mECCplexCost + "  " + mECCplexServers);
    //         //     break;
    //         // }
    //
    //     }
    // }

    private static void runJGAPGAExample() throws ClassNotFoundException, IOException, InvalidConfigurationException {
        int serverNumber = 0;
        for (int ep = 0; ep < 10; ep ++) {
            // int serverNumber = 0;
            serverNumber += 10;
            double d = 1;
            int h = 1;
            int population = 25;
            int[][] dism = GraphGenerate(serverNumber, d);
            ModelSetup(serverNumber, dism, h);
            double totalGACost = 0;
            double totalVoteCost = 0;
            for (int i = 0; i < 10; i++) {
                runJGAPGACost(population, serverNumber);
                // if (GAFitness == 0) {
                //     System.out.println("serverNumber Problem");
                //     System.out.println("false " + "serverNumber " + serverNumber + "d " + d + "h " + h);
                //     System.out.println("dism" + Arrays.deepToString(dism));
                //     // break;
                // }
                runECGreedyVoteCost();
                totalGACost += mGACost;
                totalVoteCost += mECGreedyVoteCost;
            }
            if (GAFitness == 0) break;
            double averageGACost = totalGACost / 10;
            double averageVoteCost = totalVoteCost / 10;
            System.out.println("averageGACost:  " + averageGACost);
            System.out.println("averageVoteCost:  " + averageVoteCost);

            mLines.add("==============================setting1==============================  " );
            mLines.add("ServersNumber:" + serverNumber + "  Density:" + d + "  Hops:" + h );
            mLines.add("averageGACost:  " + averageGACost);
            mLines.add("averageVoteCost:  " + averageVoteCost);
        }
        int d = 0;
        for (int ep = 0; ep < 5; ep ++) {
            // int serverNumber = 0;
            serverNumber = 150;
            // double d = 0;
            d+=1;
            int h = 1;
            int population = 25;
            int[][] dism = GraphGenerate(serverNumber, d);
            ModelSetup(serverNumber, dism, h);
            double totalGACost = 0;
            double totalVoteCost = 0;
            for (int i = 0; i < 10; i++) {
                runJGAPGACost(population, serverNumber);
                // if (GAFitness == 0) {
                //     System.out.println("density Problem");
                //     System.out.println("false " + "serverNumber " + serverNumber + "d " + d + "h " + h);
                //     System.out.println("dism" + Arrays.deepToString(dism));
                //     break;
                // }
                runECGreedyVoteCost();
                totalGACost += mGACost;
                totalVoteCost += mECGreedyVoteCost;
            }
            if (GAFitness == 0) break;

            double averageGACost = totalGACost / 10;
            double averageVoteCost = totalVoteCost / 10;
            System.out.println("averageGACost:  " + averageGACost);
            System.out.println("averageVoteCost:  " + averageVoteCost);

            mLines.add("==============================setting2==============================  " );
            mLines.add("ServersNumber:" + serverNumber + "  Density:" + d + "  Hops:" + h );
            mLines.add("averageGACost:  " + averageGACost);
            mLines.add("averageVoteCost:  " + averageVoteCost);
        }
        int h = 1;
        for (int ep = 0; ep < 5; ep ++) {
            // int serverNumber = 0;
            serverNumber = 150;
            d = 1;
            // int h = 0;
            h++;
            int population = 25;
            int[][] dism = GraphGenerate(serverNumber, d);
            ModelSetup(serverNumber, dism, h);
            double totalGACost = 0;
            double totalVoteCost = 0;
            for (int i = 0; i < 10; i++) {
                runJGAPGACost(population, serverNumber);
                if (GAFitness == 0) {
                    System.out.println("hops Problem");
                    System.out.println("false " + "serverNumber " + serverNumber + "d " + d + "h " + h);
                    System.out.println("dism" + Arrays.deepToString(dism));
                    break;
                }
                runECGreedyVoteCost();
                totalGACost += mGACost;
                totalVoteCost += mECGreedyVoteCost;
            }
            double averageGACost = totalGACost / 10;
            double averageVoteCost = totalVoteCost / 10;
            System.out.println("averageGACost:  " + averageGACost);
            System.out.println("averageVoteCost:  " + averageVoteCost);

            mLines.add("==============================setting3==============================  " );
            mLines.add("ServersNumber:" + serverNumber + "  Density:" + d + "  Hops:" + h );
            mLines.add("averageGACost:  " + averageGACost);
            mLines.add("averageVoteCost:  " + averageVoteCost);
        }
    }

    // public static void runGACost(int population, int serverNumber, ArrayList<List> cplexSolutionList) {
    //     // mGACost = mGAModel.getmGACost();
    //     long start = System.currentTimeMillis();
    //     mGAModel.runGACost(population, serverNumber, cplexSolutionList);
    //     long end = System.currentTimeMillis();
    //
    //     // System.out.println("\nGAModel");
    //     System.out.println("Time:" + (end - start) + " ms");
    //     System.out.println();
    // }

    public static void runJGAPGACost(int populationSize, int serverNumber) throws IOException, ClassNotFoundException, InvalidConfigurationException {

        // mGACost = mGAModel.getmGACost();
        long start = System.currentTimeMillis();
        mJGAPGAModel.runGACost(populationSize, serverNumber);
        long end = System.currentTimeMillis();
        mGACost = mJGAPGAModel.getCost();
        GAFitness = mJGAPGAModel.getFitness();

        // System.out.println("\nGAModel");
        // System.out.println("Time:" + (end - start) + " ms");
        // System.out.println();
    }

    private static void runECGreedyVoteCost() throws ClassNotFoundException, IOException {
        long start = System.currentTimeMillis();
        mECGreedyVoteModel.runECGreedyVoteCost();
        long end = System.currentTimeMillis();
        // double duration = (double) (Duration.between(start, end).toMillis());
        int packetsNeed = mECGreedyVoteModel.getPacketsNeed();
        mECGreedyVoteCost = mECGreedyVoteModel.getCost();
        mECGreedyVoteServers = mECGreedyVoteModel.getSelectedServers();
        mReplicaCost = mECGreedyVoteModel.getReplica();
        mReplicaServers = mECGreedyVoteModel.getSelectedReplicaServers();


        // System.out.println("\nECGreedyVoteModel");
        // System.out.println("Time:" + (end - start) + " ms");
        // System.out.println("mECGreedyVotePacketsNeeds:" + packetsNeed);
        // System.out.println("mECGreedyVoteCost:" + mECGreedyVoteCost);
        // System.out.println("mECGreedyVoteServers:" + mECGreedyVoteServers + "\n");
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

    private static void runECCplexCost() throws ClassNotFoundException, IOException {
        // System.out.println("\nECCplexModel");
        double min_cost = Double.MAX_VALUE;
        int best_pn = 0;
        List<Integer> mincostServers = new ArrayList<>();
        long start = System.currentTimeMillis();
        for (ECCplexModel it : mECCplexModels) {
            double cur_cost = 0.0;
            it.runECCplex();
            int packetsNeed = it.getPacketsNeed();
            cur_cost = it.getCost();
            List<Integer> CplexCurServers = it.getSelectedServers();
            if (packetsNeed == 1) {
                mReplicaCost = cur_cost;
                List<Integer> intermedia_replica = deepCopy(CplexCurServers);
                mReplicaServers = intermedia_replica;
                // System.out.println("Replica:" + CplexCurServers + " cost:" + cur_cost);
            }
            if (packetsNeed >= 2) {
                // System.out.println(packetsNeed + ":" + CplexCurServers + " cost:" + cur_cost);
            }
            if (packetsNeed >= 2 && cur_cost < min_cost) {
                min_cost = cur_cost;
                best_pn = packetsNeed;
                List<Integer> intermedia = deepCopy(CplexCurServers);
                mincostServers = intermedia;
            }
        }
        // System.out.println("\nBestOfCplex:" + best_pn + ":" + mincostServers + " cost:" + min_cost);
        mECCplexCost = min_cost;
        mECCplexServers = mincostServers;
        long end = System.currentTimeMillis();
        // double duration = (double) (Duration.between(start, end).toMillis());
        System.out.println("Time:" + (end - start) + " ms");
        // System.out.println("mECCplexPacketsNeeds:" + best_pn);
        System.out.println("mECCplexCost:" + mECCplexCost);
        System.out.println("mECCplexServers:" + mECCplexServers + "\n");
    }

    private static void getECCplexSolution() throws IOException, ClassNotFoundException {
        // System.out.println("\nECCplexModel");
        double min_cost = Double.MAX_VALUE;
        int best_pn = 0;
        List<Integer> mincostServers = new ArrayList<>();
        long start = System.currentTimeMillis();
        for (ECCplexModel it : mECCplexModels) {
            double cur_cost = 0.0;
            it.runECCplex();
            int packetsNeed = it.getPacketsNeed();
            cur_cost = it.getCost();
            //
            // List<Integer> M = new ArrayList<>();
            // M.add(packetsNeed);
            // List<List> solution = new ArrayList<>();
            // // solution.add(M);
            // solution.add(it.getSelectedServers());
            CplexSolutionList.add(it.getSelectedServers());
            System.out.println("packetsNeed: " + packetsNeed + ", CplexSolutionList" + it.getSelectedServers());
            // if (packetsNeed == 1) {
            //     mReplicaCost = cur_cost;
            //     List<Integer> intermedia_replica = deepCopy(CplexCurServers);
            //     mReplicaServers = intermedia_replica;
            //     // System.out.println("Replica:" + CplexCurServers + " cost:" + cur_cost);
            // }
            if (packetsNeed >= 2) {
                System.out.println(packetsNeed + ":" + CplexSolutionList + " cost:" + cur_cost);
            }
            // if (packetsNeed >= 2 && cur_cost < min_cost) {
            //     min_cost = cur_cost;
            //     best_pn = packetsNeed;
            //     List<Integer> intermedia = deepCopy(CplexCurServers);
            //     mincostServers = intermedia;
            // }
        }
        // System.out.println("\nBestOfCplex:" + best_pn + ":" + mincostServers + " cost:" + min_cost);
        mECCplexCost = min_cost;
        mECCplexServers = mincostServers;
        long end = System.currentTimeMillis();
        // double duration = (double) (Duration.between(start, end).toMillis());
        System.out.println("Time:" + (end - start) + " ms");
        // System.out.println("mECCplexPacketsNeeds:" + best_pn);
        System.out.println("mECCplexCost:" + mECCplexCost);
        System.out.println("mECCplexServers:" + mECCplexServers + "\n");
    }

    private static void getECGreedyVoteSolution() throws IOException, ClassNotFoundException {

        ArrayList<Integer> arrayList = new ArrayList<>();

        long start = System.currentTimeMillis();
        mECGreedyVoteModel.runECGreedyVoteCost();
        int mPacketNeed = mECGreedyVoteModel.getPacketsNeed();
        mECGreedyDegreeCost = mECGreedyVoteModel.getCost();
        voteSolutionList = (ArrayList<Integer>) mECGreedyVoteModel.getSelectedServers();
        mECGreedyDegreeServers = mECGreedyVoteModel.getSelectedServers();

        long end = System.currentTimeMillis();

        voteSolutionList.add(0, mPacketNeed);

        System.out.println("\nVoteModel");
        System.out.println("Time:" + (end-start) + " ms");
        System.out.println("VotePacketsNeeds:" + mPacketNeed);
        System.out.println("VoteCost:" + mECGreedyDegreeCost);
        System.out.println("VoteServers:" + mECGreedyDegreeServers + "\n");

    }



    private static void runECGreedyDegreeCost() throws ClassNotFoundException, IOException {
        long start = System.currentTimeMillis();
        mECGreedyDegreeModel.runECGreedyDegreeCost();
        long end = System.currentTimeMillis();
        // double duration = (double) (Duration.between(start, end).toMillis());
        int packetsNeed = mECGreedyDegreeModel.getPacketsNeed();
        mECGreedyDegreeCost = mECGreedyDegreeModel.getCost();
        mECGreedyDegreeServers = mECGreedyDegreeModel.getSelectedServers();
        // System.out.println("\nECGreedyDegreeModel");
        // System.out.println("Time:" + (end-start) + " ms");
        // System.out.println("mECGreedyDegreePacketsNeeds:" + packetsNeed);
        // System.out.println("mECGreedyDegreeCost:" + mECGreedyDegreeCost);
        // System.out.println("mECGreedyDegreeServers:" + mECGreedyDegreeServers + "\n");

    }

    private static void runECRandomCost() throws ClassNotFoundException, IOException {
        long start = System.currentTimeMillis();
        mECRandomModel.runECRandomCost();
        long end = System.currentTimeMillis();
        // double duration = (double) (Duration.between(start, end).toMillis());
        int packetsNeed = mECRandomModel.getPacketsNeed();
        mECRandomCost = mECRandomModel.getCost();
        mECRandomServers = mECRandomModel.getSelectedServers();

        // System.out.println("\nECRandomModel");
        // System.out.println("Time:" + (end-start) + " ms");
        // System.out.println("mECRandomPacketsNeeds:" + packetsNeed);
        // System.out.println("mECRandomCost:" + mECGreedyVoteCost);
        // System.out.println("mECRandomServers:" + mECGreedyVoteServers + "\n");
    }

    private static void runReplicaCost() throws ClassNotFoundException, IOException {
        long start = System.currentTimeMillis();
        mReplicaModel.runECRandomCost();
        long end = System.currentTimeMillis();
        // double duration = (double) (Duration.between(start, end).toMillis());
        // mReplicaCost = mReplicaModel.getReplica();
        // mReplicaServers = mReplicaModel.getSelectedReplicaServers();

        // System.out.println("\nReplicaModel");
        // System.out.println("Time:" + (end-start) + " ms");
        // System.out.println("mReplicaCost:" + mReplicaCost);
    }


}
