#include "Queue.h"

/******************************************************************************/
/*** C Headers                                                              ***/
/******************************************************************************/
#include <stdlib.h>
#include <math.h>

int MAX_NUM = 1000;         // maximum permutation steps

Queue * QueueService::queue_create()
{
	Queue * q = (Queue *) malloc(sizeof(Queue));
	q->size = 0;
	q->head = NULL;
	q->tail = NULL;
	return q;
}

void * QueueService::queue_remove(Queue * q)
{
	node * n;

	if (q->size == 0){
		return NULL;
		//free(q);//ÊÍ·ÅÄÚ´æ
	}

	n = q->head;

	if (q->size == 1)
	{
		void * cell = n->cell;
		q->size--;
		q->head = NULL;
		q->tail = NULL;
		free(n);
		return cell;
	}

	else
	{
		void * cell = n->cell;
		q->size--;
		q->head = (q->head)->next;
		free(n);
		return cell;
	}

}

void QueueService::queue_add(Queue * q, void * o)
{
	node * n = (node *) malloc(sizeof(node));


	n->cell = o;
	n->next = NULL;

	if (q->size == 0)
	{
		q->size++;
		q->head = n;
		q->tail = n;
	}

	else
	{
		q->size++;
		q->tail->next = n;
		q->tail = (q->tail)->next;
	}    
	
}

void * QueueService::queue_get(Queue * q, int index)
{
	node * n = q->head;
	void * cell = NULL;
	int i = 0;

	if (index < 0 || index >= q->size) 
		return NULL;

	for (i = 0; i < index; i++)
	{
		n = n->next;
	}

	cell = n->cell;
	return cell;
}


void * QueueService::queue_remove_index(Queue * q, int index)
{
	int i;

	if (index < 0 || index >= q->size) 
		return NULL;    

	if (q->size == 0)
		return NULL;

	if (index == 0)
	{
		return queue_remove(q);
	}

	else
	{
		node * n = q->head;
		node * p = q->head;

		void * cell = NULL;

		for (i = 0; i < index - 1; i++)
		{
			p = p->next;
		}

		n = p->next;
		cell = n->cell;
		q->size--;

		if (n == q->tail)
		{
			p->next = NULL;
			q->tail = p;
			free(n);
			return cell;
		}

		p->next = n->next;
		free(n);

		return cell;
	}
}

void QueueService::queue_permute(Queue * q)
{
	void ** cell = NULL;
	int i;
	int size = q->size;

	if (q == NULL)
		return;

	if (q->size == 0 || q->size == 1)
		return;

	cell = (void **) malloc(sizeof(void *) * q->size);

	for (i = 0; i < size && i < MAX_NUM; i++)
	{
		cell[i] = queue_remove_index(q, rand() % q->size);
		//cell[i] = queue_remove(q);
	}

	for (i = 0; i < size && i < MAX_NUM; i++)
	{
		queue_add(q, cell[i]);
	}

	free(cell);

}
