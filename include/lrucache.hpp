/**
    LRUCACHE:
    ChatGPT code modified by Roberto Grossi 2024.
**/

#include <iostream>
#include <list>
#include <unordered_map>

using namespace std;

class LRUCache {
   public:
    LRUCache(int _K) : K(_K) {}

    bool existsKey(uint64_t key) { return cacheItemsMap.find(key) != cacheItemsMap.end(); }

    double getUtility(uint64_t key) {
        auto it = cacheItemsMap.find(key);
        if (it == cacheItemsMap.end()) {
            throw runtime_error("Key not found in the LRU cache");
        } else {
            return it->second->second;
        }
    }

    void processKey(uint64_t key, double& utility) {
        auto it = cacheItemsMap.find(key);
        if (it != cacheItemsMap.end()) {
            cacheItemsList.splice(cacheItemsList.begin(), cacheItemsList, it->second);
            utility = it->second->second;
        } else {
            if (cacheItemsMap.size() == K) {
                auto last = cacheItemsList.end();
                --last;
                cacheItemsMap.erase(last->first);
                cacheItemsList.pop_back();
            }
            cacheItemsList.emplace_front(key, utility);
            cacheItemsMap[key] = cacheItemsList.begin();
        }
    }

    void printTopK() {
        for (auto it = cacheItemsList.begin(); it != cacheItemsList.end(); it++) {
            auto key = it->first;
            auto util = it->second;
            cout << key << " " << util << endl;
        }
        cout << endl;
    }

   private:
    int K;
    list<pair<uint64_t, double>> cacheItemsList;
    unordered_map<uint64_t, list<pair<uint64_t, double>>::iterator> cacheItemsMap;
};
