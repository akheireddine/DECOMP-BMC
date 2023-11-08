#pragma once

#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <vector>
#include <functional>
#include <future>
#include <assert.h>
#include <istream>
#include <algorithm>

template<typename Data>
class BlockingQueue {
private:
    std::queue<Data>                        queue;
    mutable std::mutex                      queue_mutex;
    const size_t                            queue_limit;

    bool                                    is_closed = false;

    std::condition_variable                 new_item_or_closed_event;
    std::condition_variable                 item_removed_event;

#ifndef NDEBUG
    size_t                                  pushes_in_progress = 0;
#endif

public:
    BlockingQueue(size_t size_limit=0) : queue_limit(size_limit)
    {}

    void push(const Data& data)
    {
        // boost::mutex::scoped_lock lock(queue_mutex);
        std::unique_lock<std::mutex> unique_queue(queue_mutex);
#ifndef NDEBUG
        ++pushes_in_progress;
#endif
        if (queue_limit > 0) {
            while (queue.size() >= queue_limit) {
                item_removed_event.wait(unique_queue);
            }
        }
        assert (!is_closed);
        queue.push(data);
#ifndef NDEBUG
        --pushes_in_progress;
#endif
        // lock.unlock();
        unique_queue.unlock();
        new_item_or_closed_event.notify_one();
    }

    bool try_push(const Data& data)
    {
        // boost::mutex::scoped_lock lock(queue_mutex);
        std::unique_lock<std::mutex> unique_queue(queue_mutex);
        if (queue_limit > 0) {
            if (queue.size() >= queue_limit) {
                return false;
            }
        }
        assert (!is_closed);
        queue.push(data);
        // lock.unlock();
        unique_queue.unlock();

        new_item_or_closed_event.notify_one();
        return true;
    }

    void close()
    {
        // boost::mutex::scoped_lock lock(queue_mutex);
        std::unique_lock<std::mutex> unique_queue(queue_mutex);
        assert (!is_closed);
#ifndef NDEBUG
        assert (pushes_in_progress == 0);
#endif
        is_closed = true;
        // lock.unlock();
        unique_queue.unlock();

        new_item_or_closed_event.notify_all();
    }

    bool pop(Data &popped_value)
    {
        // boost::mutex::scoped_lock lock(queue_mutex);
        std::unique_lock<std::mutex> unique_queue(queue_mutex);
        while (queue.empty()) {
            if (is_closed) {
                return false;
            }
            new_item_or_closed_event.wait(unique_queue);
        }

        popped_value = queue.front();
        queue.pop();
        item_removed_event.notify_one();
        return true;
    }

    bool try_pop(Data &popped_value)
    {
        // boost::mutex::scoped_lock lock(queue_mutex);
        std::unique_lock<std::mutex> unique_queue(queue_mutex);
        if (queue.empty()) {
            return false;
        }

        popped_value = queue.front();
        queue.pop();
        item_removed_event.notify_one();
        return true;
    }

    bool empty() const
    {
        // boost::mutex::scoped_lock lock(queue_mutex);
        std::unique_lock<std::mutex> unique_queue(queue_mutex);
        return queue.empty();
    }

    bool closed() const
    {
        // boost::mutex::scoped_lock lock(queue_mutex);
        std::unique_lock<std::mutex> unique_queue(queue_mutex);
        return is_closed;
    }

    size_t limit() const
    {
        return queue_limit;
    }

    size_t size() const
    {
        // boost::mutex::scoped_lock lock(queue_mutex);
        std::unique_lock<std::mutex> unique_queue(queue_mutex);
        return queue.size();
    }

};



class thread_pool
{
public:
	thread_pool(unsigned int threads = std::thread::hardware_concurrency())
	: m_queues(threads), m_count(threads)
	{
      m_index = 0;
		assert(threads != 0);
		auto worker = [&](unsigned int i)
		{
			while(true)
			{
				Proc f;
				for(unsigned int n = 0; n < m_count; n++)
					if(m_queues[(i + n) % m_count].try_pop(f)) break;
				if(!f && !m_queues[i].pop(f)) break;
				f();
			}
		};
		for(unsigned int i = 0; i < threads; ++i)
			m_threads.emplace_back(worker, i);
	}
 
	~thread_pool() noexcept
	{
		for(auto& queue : m_queues)
			queue.close();
			// queue.done();
		for(auto& thread : m_threads)
			thread.join();
	}
 
	template<typename F, typename... Args>
	void enqueue_work(F&& f, Args&&... args)
	{
		auto work = [f,args...]() { f(args...); };
		unsigned int i = m_index++;

		// for(unsigned int n = 0; n < m_count; ++n)
            // if(m_queues[n].empty() &&  m_queues[n].try_push(work) ) return;
		for(unsigned int n = 0; n < m_count * K; ++n)
			if(m_queues[(i + n) % m_count].try_push(work)) return;
		m_queues[i % m_count].push(work);
	}
 
	template<typename F, typename... Args>
	auto enqueue_task(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>
	{
		using return_type = typename std::result_of<F(Args...)>::type;
		auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
		std::future<return_type> res = task->get_future();
 
		auto work = [task](){ (*task)(); };
		unsigned int i = m_index++;
		for(unsigned int n = 0; n < m_count * K; ++n)
			if(m_queues[(i + n) % m_count].try_push(work)) return res;
		m_queues[i % m_count].push(work);
 
		return res;
	}

    int get_remains_thread()
    {
        return m_queues.size();
    }
 
// protected:
	using Proc = std::function<void(void)>;
	using Queues = std::vector<BlockingQueue<Proc>>;
	Queues m_queues;
 
	using Threads = std::vector<std::thread>;
	Threads m_threads;
 
	const unsigned int m_count;
	std::atomic_uint m_index;
 
	static const unsigned int K = 3;
};