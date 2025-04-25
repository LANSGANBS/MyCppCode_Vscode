#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>

std::random_device rd;
std::mt19937 gen(rd());

// [min, max] 上的随机数生成器
int64_t randomInt(int64_t min, int64_t max) {
  std::uniform_int_distribution<int64_t> dis(min, max);
  return dis(gen);
}

std::pair<std::string, int64_t> run(const std::string& path,
                                    const std::string& input) {
  const std::string tempInputFile = "input.txt";
  const std::string tempOutputFile = "output.txt";

  // 写入 input.txt（不存在时创建，存在时覆盖）
  std::ofstream inputFile(tempInputFile);
  inputFile << input;
  inputFile.close();

  // 执行命令，将结果重定向到 output.txt（同样会创建或覆盖文件）23
  std::string command = path + " < " + tempInputFile + " > " + tempOutputFile;
  auto start = std::chrono::high_resolution_clock::now();
  system(command.c_str());
  auto end = std::chrono::high_resolution_clock::now();
  int64_t duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();

  // 读取 output.txt
  std::ifstream outputFile(tempOutputFile);
  std::string output((std::istreambuf_iterator<char>(outputFile)),
                     std::istreambuf_iterator<char>());
  outputFile.close();

  // 不再删除 input.txt
  // remove(tempInputFile.c_str());
  remove(tempOutputFile.c_str());

  return {output, duration};
}

// 数据生成函数
std::string gen_input() {
  // 每次生成一个测试文件的内容, 请自行设计
  std::ostringstream oss;
  //   int n = randomInt(1, 200);
  //   oss << n << "\n";
  //   for (int i = 0; i < n; ++i) {
  //     int x = randomInt(-200, 200);
  //     oss << x << " \n"[i == n - 1];
  //   }
  //   int q = randomInt(1, 200);
  //   oss << q << "\n";
  //   for (int i = 0; i < q; ++i) {
  //     int l = randomInt(1, n);
  //     int r = randomInt(l, n);
  //     oss << l << " " << r << "\n";
  //   }
  // ---
  int a = randomInt(1, 200), b = randomInt(1, 200);
  oss << a << ' ' << b;
  return oss.str();
}

int main() {
  std::string path1 = "1.exe";
  std::string path2 = "2.exe";

  // 对拍次数
  int64_t testCases = 10000;

  for (int64_t i = 1; i <= testCases; ++i) {
    std::string input = gen_input();
    auto res1 = run(path1, input);
    auto res2 = run(path2, input);
    std::string output1 = res1.first;
    std::string output2 = res2.first;
    int64_t duration1 = res1.second;
    int64_t duration2 = res2.second;

    std::cout << "Test case " << i << '\n'
              << "Input: \n"
              << input << '\n'
              << "Output 1: \n"
              << output1 << '\n'
              << "Output 2: \n"
              << output2 << '\n'
              << "Time 1: " << duration1 << " ms" << '\n'
              << "Time 2: " << duration2 << " ms" << '\n';

    if (output1 != output2) {
      std::cout << "Test case " << i << " failed!" << '\n';
      return 1;
    }

    std::cout << "Test case " << i << " passed!" << '\n';
  }

  std::cout << "All test cases passed!" << '\n';
  return 0;
}
