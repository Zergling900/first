#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>   // isalpha

using namespace std;

struct Coord {
    double x, y, z;
};

bool is_comment_or_empty(const string &line) {
    // 去掉前导空格
    size_t i = 0;
    while (i < line.size() && isspace(static_cast<unsigned char>(line[i]))) i++;

    if (i == line.size()) return true;                // 空行
    if (line[i] == '#' || line[i] == '/') return true; // 注释行（# 或 // 开头）
    return false;
}

int main() {
    ifstream fin("input.txt");
    if (!fin) {
        cerr << "Cannot open file\n";
        return 1;
    }

    string key;
    vector<Coord> coords;

    // 先正常按 key-value 方式读，找到 "coordinates"
    while (fin >> key) {
        if (key == "coordinates") {
            // 读掉本行剩余部分的换行符
            string dummy;
            getline(fin, dummy);

            // 接下来逐行读坐标
            string line;
            while (getline(fin, line)) {
                if (is_comment_or_empty(line)) break;  // 碰到空行或注释就停止这个 coordinates 区块

                stringstream ss(line);
                Coord c;
                if (ss >> c.x >> c.y >> c.z) {
                    coords.push_back(c);
                } else {
                    // 如果这一行不是三个数，说明到了下一段内容，回退流位置或者直接 break
                    break;
                }
            }
        }
        // 其他 key 就按你原来的方式去处理，比如 energy, Box_Lx 等
    }

    // 检查一下是否读到了全部坐标
    cout << "Total coordinates: " << coords.size() << "\n";
    for (size_t i = 0; i < coords.size(); ++i) {
        cout << i << ": " << coords[i].x << " " << coords[i].y << " " << coords[i].z << "\n";
    }

    return 0;
}
