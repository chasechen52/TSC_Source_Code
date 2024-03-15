package tool;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;

public class DistanceMatrixManager {
    private static DistanceMatrixManager instance;
    private static HashMap<String, int[][]> distanceMatrices = new HashMap<>();
    private static String filePath;

    private DistanceMatrixManager(String filePath) {
        this.filePath = filePath;
        loadMatricesFromFile();
    }

    public static DistanceMatrixManager getInstance(String filePath) {
        if (instance == null) {
            instance = new DistanceMatrixManager(filePath);
        }
        return instance;
    }

    /**
     * 存储距离矩阵
     * @param key 矩阵的键值,用于后续获取
     * @param matrix 距离矩阵
     */
    public void storeDistanceMatrix(String key, int[][] matrix) {
        distanceMatrices.put(key, matrix);
        saveMatricesToFile();
    }

    /**
     * 获取距离矩阵
     * @param key 矩阵的键值
     * @return 距离矩阵,如果不存在则返回null
     */
    public int[][] getDistanceMatrix(String key) {
        // System.out.println(Arrays.deepToString(distanceMatrices.get(key)));
        return distanceMatrices.get(key);
    }

    public void loadMatricesFromFile() {
        try {
            FileInputStream fis = new FileInputStream(filePath);
            ObjectInputStream ois = new ObjectInputStream(fis);
            distanceMatrices = (HashMap<String, int[][]>) ois.readObject();
            ois.close();
            fis.close();
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
        }
    }

    private void saveMatricesToFile() {
        try {
            FileOutputStream fos = new FileOutputStream(filePath);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(distanceMatrices);
            oos.flush();
            oos.close();
            fos.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}