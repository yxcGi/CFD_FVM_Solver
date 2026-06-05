#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <iostream>
#include <thread>
#include <future>
#include <queue>
#include <atomic>
#include <functional>
#include <unordered_map>
#include <chrono>


static constexpr int DEFAULT_INIT_THREAD_NUM = 4;      // 默认初始线程数量
static constexpr int MAX_TASK_QUEUE_NUM = 100;         // 最大任务队列数量


class ThreadPool
{
public:
    ThreadPool() = default;

    explicit ThreadPool(std::size_t threadNum,
                        std::size_t maxTaskQueNum = MAX_TASK_QUEUE_NUM)
    {
        start(threadNum, maxTaskQueNum);
    }

    ~ThreadPool()
    {
        shutdown();
    }

    static ThreadPool& getInstance()
    {
        static ThreadPool instance;
        return instance;
    }

    void start(std::size_t threadNum = DEFAULT_INIT_THREAD_NUM,
               std::size_t maxTaskQueNum = MAX_TASK_QUEUE_NUM)
    {
        if (threadNum == 0)
        {
            threadNum = 1;
        }

        std::unique_lock<std::mutex> lock(mtx_);

        // 防止重复 start() 反复创建线程。
        if (running_)
        {
            return;
        }

        running_ = true;
        stopRequested_ = false;
        maxTaskQueNum_ = std::max<std::size_t>(1, maxTaskQueNum);

        workers_.reserve(threadNum);

        for (std::size_t i = 0; i < threadNum; ++i)
        {
            workers_.emplace_back(&ThreadPool::workerLoop, this);
        }
    }

    void shutdown() noexcept
    {
        {
            std::unique_lock<std::mutex> lock(mtx_);

            if (!running_)
            {
                return;
            }

            stopRequested_ = true;
        }

        notEmpty_.notify_all();

        for (std::thread& worker : workers_)
        {
            if (worker.joinable())
            {
                worker.join();
            }
        }

        {
            std::unique_lock<std::mutex> lock(mtx_);
            workers_.clear();
            running_ = false;
            stopRequested_ = false;
        }
    }

    template<typename Func, typename... Args>
    auto submitTask(Func&& func, Args&&... args)
    {
        using ReturnType = std::invoke_result_t<Func, Args...>;

        auto boundTask = std::bind(
            std::forward<Func>(func),
            std::forward<Args>(args)...
        );

        auto packagedTask =
            std::make_shared<std::packaged_task<ReturnType()>>(
                std::move(boundTask)
            );

        std::future<ReturnType> result = packagedTask->get_future();

        {
            std::unique_lock<std::mutex> lock(mtx_);

            if (!running_ || stopRequested_)
            {
                throw std::runtime_error("ThreadPool is not running.");
            }

            notFull_.wait(lock, [this]() {
                return stopRequested_ || taskQue_.size() < maxTaskQueNum_;
                          });

            if (stopRequested_)
            {
                throw std::runtime_error("ThreadPool is stopped.");
            }

            taskQue_.emplace([packagedTask]() {
                (*packagedTask)();
                             });
        }

        notEmpty_.notify_one();
        return result;
    }

private:
    void workerLoop()
    {
        while (true)
        {
            std::function<void()> task;

            {
                std::unique_lock<std::mutex> lock(mtx_);

                notEmpty_.wait(lock, [this]() {
                    return stopRequested_ || !taskQue_.empty();
                               });

                if (stopRequested_ && taskQue_.empty())
                {
                    return;
                }

                task = std::move(taskQue_.front());
                taskQue_.pop();
            }

            notFull_.notify_one();
            task();
        }
    }

private:
    mutable std::mutex mtx_;
    std::condition_variable notFull_;
    std::condition_variable notEmpty_;

    std::queue<std::function<void()>> taskQue_;
    std::vector<std::thread> workers_;

    std::size_t maxTaskQueNum_{ MAX_TASK_QUEUE_NUM };
    bool running_{ false };
    bool stopRequested_{ false };
};
#endif // THREADPOOL_H_