package bam

import "container/list"

type blockCache interface {
	Get(key int64) ([]byte, bool)
	Set(key int64, value []byte)
}

////
// uses maps for storing small data sets

type mapCache map[int64][]byte

func newMapCache(n int) blockCache {
	return mapCache(make(map[int64][]byte, n))
}

func (x mapCache) Get(key int64) ([]byte, bool) {
	r, b := x[key]
	return r, b
}

func (x mapCache) Set(key int64, value []byte) {
	x[key] = value
}

////
// uses a more advanced cache for storing large data sets
// S4-LRU described in http://www.cs.cornell.edu/~qhuang/papers/sosp_fbanalysis.pdf

type cacheItem struct {
	qid   int
	key   int64
	value []byte
}

type blockLRUCache struct {
	cap    int
	data   map[int64]*list.Element
	queues []*list.List
}

func newLRUCache(capacity int) blockCache {
	return &blockLRUCache{
		cap:    (capacity + 3) / 4,
		data:   make(map[int64]*list.Element),
		queues: []*list.List{list.New(), list.New(), list.New(), list.New()},
	}
}

func (c *blockLRUCache) Get(key int64) ([]byte, bool) {
	v, ok := c.data[key]
	if !ok {
		return nil, false
	}

	item := v.Value.(*cacheItem)

	if item.qid == 3 {
		// can't bump up a level
		c.queues[3].MoveToFront(v)
		return item.value, true
	}

	if c.queues[item.qid+1].Len() < c.cap {
		// plenty of room to bump up
		c.queues[item.qid].Remove(v)
		item.qid++
		c.data[key] = c.queues[item.qid].PushFront(item)
		return item.value, true
	}

	// full queues, swap items in place
	w := c.queues[item.qid+1].Back()
	other := w.Value.(*cacheItem)

	// keep qid on each
	other.key, item.key = item.key, other.key
	other.value, item.value = item.value, other.value

	c.data[item.key] = v
	c.data[other.key] = w
	c.queues[item.qid].MoveToFront(v)
	c.queues[other.qid].MoveToFront(w)

	return other.value, true
}

func (c *blockLRUCache) Set(key int64, value []byte) {
	if c.queues[0].Len() < c.cap {
		c.data[key] = c.queues[0].PushFront(&cacheItem{0, key, value})
		return
	}

	// reuse the tail item
	e := c.queues[0].Back()
	item := e.Value.(*cacheItem)

	delete(c.data, item.key)
	item.key = key
	item.value = value
	c.data[key] = e
	c.queues[0].MoveToFront(e)
}
