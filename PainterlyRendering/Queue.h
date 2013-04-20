#ifndef _QUEUE_H_
#define _QUEUE_H_

#pragma once

typedef struct node
{
	void * cell;
	struct node * next;

} node;

typedef struct Queue
{
	int size;
	node * head;
	node * tail;
} Queue;

class QueueService
{
public:
	QueueService(void){}
	~QueueService(void){}

	static Queue * queue_create();

	static void * queue_remove(Queue * q);

	static void queue_add(Queue * q, void * o);

	static void * queue_get(Queue * q, int index);

	static void queue_permute(Queue * q);

	static void * queue_remove_index(Queue * q, int index);
};

#endif
