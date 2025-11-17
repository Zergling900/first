#include <string>
#include <sstream>
#include <cctype>   // isalpha
#include <algorithm>
#include <fstream>


using namespace std;

// 清理一行：去掉注释、去掉空格
std::string clean_line(std::string line)
{
    // 去掉行尾整段注释：// 注释
    {
        size_t pos = line.find("//");
        if (pos != std::string::npos)
            line = line.substr(0, pos);
    }

    // 去掉行尾 # 注释（但不删 key 中间的 # 例如 a#b）
    {
        size_t pos = line.find("#");
        if (pos != std::string::npos)
            line = line.substr(0, pos);
    }

    // 去掉单行格式的 /* ... */
    {
        size_t pos = line.find("/*");
        if (pos != std::string::npos)
        {
            size_t end = line.find("*/", pos + 2);
            if (end != std::string::npos)
                line.erase(pos, end - pos + 2);
            else
                line = line.substr(0, pos);
        }
    }

    // 去掉前后空白字符
    auto not_space = [](int ch){ return !std::isspace(ch); };

    // 去掉前空白
    line.erase(line.begin(),
               std::find_if(line.begin(), line.end(), not_space));

    // 如果删掉后为空，则直接返回空行
    if (line.empty()) return "";

    // 去掉后空白
    line.erase(std::find_if(line.rbegin(), line.rend(), not_space).base(),
               line.end());

    // 检查行首是否为注释符号 (# 或 // 或 /*)
    if (line.rfind("//", 0) == 0) return "";
    if (line.rfind("#",  0) == 0) return "";
    if (line.rfind("/*", 0) == 0) return "";

    return line;
}

bool read_clean_line(std::ifstream &fin, std::string &out)
{
    std::string line;

    while (std::getline(fin, line))
    {
        line = clean_line(line);

        if (line.empty())  // 空行 or 注释行
            continue;

        out = line;
        return true; // 找到有效内容
    }

    return false;  // EOF
}
