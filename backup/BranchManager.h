#ifndef BRANCHMANAGER_H
#define BRANCHMANAGER_H

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>

#include "TBranch.h"
#include "TTree.h"

class BranchBase {
public:
    explicit BranchBase(const char *name) : branchName(name) {}
    virtual ~BranchBase() = default;

    virtual void attach(TTree *tree_) {
        tree = tree_;
        branch = tree ? tree->GetBranch(branchName.c_str()) : nullptr;
        const bool wasActive = active; // remember whether this branch had been turned on
        lastEntry = -1;
        needsRebind = branch != nullptr;
        active = wasActive && branch;
        if (tree && branch) {
            tree->SetBranchStatus(branchName.c_str(), wasActive ? 1 : 0);
            if (active)
                bindAddress();
        }
    }

    static void recordActiveBranch(const std::string &name) {
        activeBranchNames.insert(name);
    }

    static const std::unordered_set<std::string> &GetActiveBranches() {
        return activeBranchNames;
    }

    static void ClearActiveBranches() { activeBranchNames.clear(); }

    virtual void reset() { lastEntry = -1; }

    void bindEntrySource(const Long64_t *ptr) { entrySource = ptr; }

protected:
    virtual void bindAddress() const {}

    Long64_t resolveEntry(Long64_t entry) const {
        if (entry >= 0)
            return entry;
        return entrySource ? *entrySource : -1;
    }

    std::string branchName;
    TTree *tree = nullptr;
    TBranch *branch = nullptr;
    mutable bool active = false;
    mutable Long64_t lastEntry = -1;
    const Long64_t *entrySource = nullptr;
    mutable bool needsRebind = false;

    inline static std::unordered_set<std::string> activeBranchNames{};
};

template <typename T>
class BranchScalar : public BranchBase {
public:
    explicit BranchScalar(const char *name) : BranchBase(name) {}

    void bindAddress() const override {
        if (branch)
        {
            branch->SetAddress(&value);
            needsRebind = false;
        }
    }

    const T &get(Long64_t entry = -1) const {
        ensure(entry);
        return value;
    }

    operator T() const { return get(); }

    T &mutableValue() const { return value; }

    bool valid() const { return branch != nullptr; }

    void ensure(Long64_t entry = -1) const {
        if (!branch)
            return;
        if (needsRebind)
            bindAddress();
        const Long64_t target = resolveEntry(entry);
        if (target < 0)
            return;
        if (!active) {
            std::cout << "[BranchManager] activating branch '" << branchName
                      << "'" << std::endl;
            tree->SetBranchStatus(branchName.c_str(), 1);
            branch->SetAddress(&value);
            active = true;
            needsRebind = false;
            recordActiveBranch(branchName);
        }
        if (lastEntry == target)
            return;
        branch->GetEntry(target);
        lastEntry = target;
    }

private:
    mutable T value{};
};

template <typename T, typename CountT>
class BranchVector : public BranchBase {
public:
    BranchVector(const char *name, BranchScalar<CountT> &count)
        : BranchBase(name), countBranch(count) {}

    void bindAddress() const override {
        if (!branch)
            return;
        if (capacity > 0 && !buffer.empty())
            branch->SetAddress(buffer.data());
        else
            branch->SetAddress(&zeroValue);
        needsRebind = false;
    }

    const std::vector<T> &values(Long64_t entry = -1) const {
        ensure(entry);
        return buffer;
    }

    std::vector<T> &mutableValues(Long64_t entry = -1) const {
        ensure(entry);
        return buffer;
    }

    const T &operator[](std::size_t idx) const {
        ensure();
        return buffer.at(idx);
    }

    std::size_t size(Long64_t entry = -1) const {
        ensure(entry);
        return buffer.size();
    }

    bool empty(Long64_t entry = -1) const {
        return size(entry) == 0;
    }

    bool valid() const { return branch != nullptr; }

    void ensure(Long64_t entry = -1) const {
        if (!branch)
            return;
        if (needsRebind)
            bindAddress();
        if (!active) {
            std::cout << "[BranchManager] activating branch '" << branchName
                      << "'" << std::endl;
            tree->SetBranchStatus(branchName.c_str(), 1);
            active = true;
            capacity = 0;
            recordActiveBranch(branchName);
        }

        const Long64_t target = resolveEntry(entry);
        if (target < 0)
            return;

        countBranch.ensure(target);
        auto rawSize = countBranch.mutableValue();
        if (rawSize < 0)
            rawSize = 0;
        std::size_t needed = static_cast<std::size_t>(rawSize);

        if (needed == 0)
        {
            buffer.clear();
            branch->SetAddress(&zeroValue);
            capacity = 0;
            if (lastEntry != target)
            {
                branch->GetEntry(target);
                lastEntry = target;
            }
            return;
        }

        if (buffer.capacity() < needed)
            buffer.reserve(needed);
        buffer.resize(needed);
        branch->SetAddress(buffer.data());
        capacity = buffer.capacity();

        if (lastEntry == target)
            return;
        branch->GetEntry(target);
        lastEntry = target;
    }

private:
    BranchScalar<CountT> &countBranch;
    mutable std::vector<T> buffer;
    mutable std::size_t capacity = 0;
    mutable T zeroValue{};
};

template <typename CountT>
class BranchVector<bool, CountT> : public BranchBase {
public:
    BranchVector(const char *name, BranchScalar<CountT> &count)
        : BranchBase(name), countBranch(count) {}

    void bindAddress() const override {
        if (!branch)
            return;
        if (capacity > 0 && !storage.empty())
            branch->SetAddress(storage.data());
        else
            branch->SetAddress(&zeroValue);
        needsRebind = false;
    }

    bool operator[](std::size_t idx) const {
        ensure();
        return static_cast<bool>(storage.at(idx));
    }

    std::size_t size(Long64_t entry = -1) const {
        ensure(entry);
        return storage.size();
    }

    bool empty(Long64_t entry = -1) const {
        return size(entry) == 0;
    }

    bool valid() const { return branch != nullptr; }

    void ensure(Long64_t entry = -1) const {
        if (!branch)
            return;
        if (needsRebind)
            bindAddress();
        if (!active) {
            tree->SetBranchStatus(branchName.c_str(), 1);
            std::cout << "[BranchManager] activating branch '" << branchName << "'" << std::endl;
            active = true;
            capacity = 0;
            recordActiveBranch(branchName);
        }

        const Long64_t target = resolveEntry(entry);
        if (target < 0)
            return;

        countBranch.ensure(target);
        auto rawSize = countBranch.mutableValue();
        if (rawSize < 0)
            rawSize = 0;
        std::size_t needed = static_cast<std::size_t>(rawSize);

        if (needed == 0)
        {
            storage.clear();
            branch->SetAddress(&zeroValue);
            capacity = 0;
            if (lastEntry != target)
            {
                branch->GetEntry(target);
                lastEntry = target;
            }
            return;
        }

        if (storage.capacity() < needed)
            storage.reserve(needed);
        storage.resize(needed);
        branch->SetAddress(storage.data());
        capacity = storage.capacity();

        if (lastEntry == target)
            return;
        branch->GetEntry(target);
        lastEntry = target;
    }

private:
    BranchScalar<CountT> &countBranch;
    mutable std::vector<unsigned char> storage;
    mutable std::size_t capacity = 0;
    mutable unsigned char zeroValue = 0;
};

class BranchManager {
public:
    void clear() {
        branches.clear();
        tree = nullptr;
    }

    void bindEntrySource(const Long64_t *ptr) {
        entrySource = ptr;
        for (auto *branch : branches)
            branch->bindEntrySource(ptr);
    }

    void attachTree(TTree *t) {
        tree = t;
        for (auto *branch : branches) {
            branch->attach(tree);
            branch->bindEntrySource(entrySource);
        }
    }

    template <typename Branch>
    void registerBranch(Branch &branch) {
        branch.reset();
        branch.bindEntrySource(entrySource);
        branches.push_back(&branch);
        if (tree)
            branch.attach(tree);
    }

    template <typename T>
    void registerScalar(BranchScalar<T> &branch) { registerBranch(branch); }

    template <typename T, typename CountT>
    void registerVector(BranchVector<T, CountT> &branch) { registerBranch(branch); }

    template <typename CountT>
    void registerVector(BranchVector<bool, CountT> &branch) { registerBranch(branch); }

    template <typename Branch>
    void reset(Branch &branch) { branch.reset(); }

    void resetAll() {
        for (auto *branch : branches)
            branch->reset();
    }

private:
    TTree *tree = nullptr;
    const Long64_t *entrySource = nullptr;
    std::vector<BranchBase *> branches;
};

#endif // BRANCHMANAGER_H
