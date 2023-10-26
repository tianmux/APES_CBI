#include <jsoncpp/json/json.h>
#include <fstream>
#include <iostream>
int main() {
    Json::Value root;
    Json::CharReaderBuilder builder;
    std::ifstream config_doc("config.json", std::ifstream::binary);
    std::string errs;
    bool ok = Json::parseFromStream(builder, config_doc, &root, &errs);
    if(!ok) {
        // Handle error
    }
    for (int i=0; i<root.size(); i++) {
        std::cout << root.getMemberNames()[i]<<":";
        for (int j=0; j<root[root.getMemberNames()[i]].size(); j++) {
            std::cout << root[root.getMemberNames()[i]][j] << ",";
        }
        std::cout << std::endl;
    }
    return 0;
}