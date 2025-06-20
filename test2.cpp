#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

// 计算给定顺序的磁道移动距离，同时记录经过的磁道序列。
// currentTrack 表示起始磁道号，tracks 表示接下来的磁道访问顺序。
pair<int, vector<int>> calculateMovement(const vector<int>& tracks,
                                         int currentTrack) {
  int totalMovement = 0;
  vector<int> visitedTracks;
  visitedTracks.push_back(currentTrack);
  for (int track : tracks) {
    totalMovement += abs(currentTrack - track);
    currentTrack = track;
    visitedTracks.push_back(currentTrack);
  }
  return {totalMovement, visitedTracks};
}

// FIFO算法：按照给定顺序访问磁道
pair<int, vector<int>> fifo(const vector<int>& tracks, int currentTrack) {
  return calculateMovement(tracks, currentTrack);
}

// SSTF算法：每次选择距离当前磁道最近的磁道进行访问
pair<int, vector<int>> sstf(const vector<int>& tracks, int currentTrack) {
  vector<int> remaining = tracks;
  vector<int> visitedTracks;
  visitedTracks.push_back(currentTrack);
  int totalMovement = 0;
  while (!remaining.empty()) {
    int closestIndex = 0;
    int closestDistance = abs(currentTrack - remaining[0]);
    for (int i = 0; i < (int)remaining.size(); i++) {
      int distance = abs(currentTrack - remaining[i]);
      if (distance < closestDistance) {
        closestDistance = distance;
        closestIndex = i;
      }
    }
    totalMovement += closestDistance;
    currentTrack = remaining[closestIndex];
    visitedTracks.push_back(currentTrack);
    remaining.erase(remaining.begin() + closestIndex);
  }
  return {totalMovement, visitedTracks};
}

// SCAN算法：先向一个方向扫描，再折返扫描另一侧（先访问大于等于起始磁道的，再逆向访问低于起始磁道的）
pair<int, vector<int>> scan(const vector<int>& tracks, int currentTrack) {
  vector<int> sorted = tracks;
  sort(sorted.begin(), sorted.end());
  vector<int> left, right;
  for (int t : sorted) {
    if (t < currentTrack)
      left.push_back(t);
    else
      right.push_back(t);
  }
  // 先访问右侧，再逆转左侧后访问
  vector<int> visitOrder = right;
  reverse(left.begin(), left.end());
  visitOrder.insert(visitOrder.end(), left.begin(), left.end());
  return calculateMovement(visitOrder, currentTrack);
}

// C-SCAN算法：单向扫描，扫描到端后返回另一端接着扫描
pair<int, vector<int>> cscan(const vector<int>& tracks, int currentTrack) {
  vector<int> sorted = tracks;
  sort(sorted.begin(), sorted.end());
  vector<int> left, right;
  for (int t : sorted) {
    if (t < currentTrack)
      left.push_back(t);
    else
      right.push_back(t);
  }
  // 先访问右侧，再访问左侧（模拟回到起点的动作忽略回退距离）
  vector<int> visitOrder = right;
  visitOrder.insert(visitOrder.end(), left.begin(), left.end());
  return calculateMovement(visitOrder, currentTrack);
}

// FSCAN算法：本题中与SCAN的实现一致
pair<int, vector<int>> fscan(const vector<int>& tracks, int currentTrack) {
  return scan(tracks, currentTrack);
}

int main() {
  // 使用随机数生成器生成磁道号。范围为 [0, 199]
  mt19937 rng((unsigned)time(nullptr));
  uniform_int_distribution<int> dist(0, 199);

  int numTracks;
  cout << "请输入要处理的磁道数: ";
  cin >> numTracks;
  vector<int> tracks(numTracks);
  for (int i = 0; i < numTracks; i++) {
    tracks[i] = dist(rng);
  }
  cout << "随机生成的磁道号为: ";
  for (int t : tracks) {
    cout << t << " ";
  }
  cout << "\n";

  while (true) {
    cout << "\n请选择算法:\n";
    cout << "1. FIFO\n";
    cout << "2. SSTF\n";
    cout << "3. SCAN\n";
    cout << "4. C-SCAN\n";
    cout << "5. FSCAN\n";
    cout << "请输入选择 (1-5): ";
    int choice;
    cin >> choice;

    cout << "请输入当前磁道号: ";
    int currentTrack;
    cin >> currentTrack;

    vector<int> sortedTracks = tracks;  // 拷贝原始磁道号
    // 针对不同选择进行预处理排序
    if (choice == 1) {
      // FIFO：按照原始顺序
      // sortedTracks 已经赋值，无需排序
    } else if (choice == 2) {
      // SSTF：按距离当前磁道排序，距离较近的靠前
      sort(sortedTracks.begin(), sortedTracks.end(),
           [currentTrack](int a, int b) {
             return abs(a - currentTrack) < abs(b - currentTrack);
           });
    } else if (choice >= 3 && choice <= 5) {
      // SCAN, C-SCAN, FSCAN：升序排序
      sort(sortedTracks.begin(), sortedTracks.end());
    } else {
      cout << "无效选择，请重新输入。\n";
      continue;
    }

    cout << "排序后的磁道分布: ";
    for (int t : sortedTracks) {
      cout << t << " ";
    }
    cout << "\n";

    pair<int, vector<int>> res;
    switch (choice) {
      case 1:
        res = fifo(sortedTracks, currentTrack);
        break;
      case 2:
        res = sstf(sortedTracks, currentTrack);
        break;
      case 3:
        res = scan(sortedTracks, currentTrack);
        break;
      case 4:
        res = cscan(sortedTracks, currentTrack);
        break;
      case 5:
        res = fscan(sortedTracks, currentTrack);
        break;
      default:
        cout << "无效选择，请重新输入。\n";
        continue;
    }

    int totalMovement = res.first;
    // 访问的磁道数包含起始点
    int visitedCount = res.second.size() - 1;  // 中间移动次数
    double avgMovement =
        visitedCount > 0 ? (double)totalMovement / visitedCount : 0.0;

    cout << "移动的总磁道数: " << totalMovement << "\n";
    cout << "经过的总磁道数: " << visitedCount << "\n";
    cout << "移动的平均磁道数: " << fixed << setprecision(2) << avgMovement
         << "\n";

    cout << "\n是否再次选择算法? (y/n): ";
    string answer;
    cin >> answer;
    if (answer != "y" && answer != "Y") break;
  }
  return 0;
}