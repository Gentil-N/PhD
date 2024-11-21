#ifndef __CUTILS_BASE_H__
#define __CUTILS_BASE_H__

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*
*   MISCELLANEOUS
*/
#define for_loop(index_name, limit) for (size_t index_name = 0; index_name < limit; ++index_name)
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))
#define clamp(x, low, up) min((up), max((x), (low)))

/*
*   STRING
*/
#define filename(path) strrchr(path, '/') ? strrchr(path, '/') + 1 : path

/*
*   LOGGER
*/

#define LOG_DEFINE(suffix, type)                                                                                    \
typedef void (*suffix##_##type##_callback_fn)(int level, int code, const char *file, int line, const char *msg);    \
int suffix##_set_##type##_callback(suffix##_##type##_callback_fn fn);                                               \
int _##suffix##_throw_##type(int level, int code, const char *file, int line, const char *fmt, ...);                \

#define LOG_DECLARE(suffix, type)                                                                       \
static int suffix##_##type##_INIT = 0;                                                                  \
static pthread_mutex_t suffix##_##type##_MUTEX = {0};                                                   \
static suffix##_##type##_callback_fn suffix##_##type##_CALLBACK = NULL;                                 \
int suffix##_set_##type##_callback(suffix##_##type##_callback_fn fn)                                    \
{                                                                                                       \
    if (fn == NULL && suffix##_##type##_INIT == 1)                                                      \
    {                                                                                                   \
        suffix##_##type##_INIT = 0;                                                                     \
        suffix##_##type##_CALLBACK = NULL;                                                              \
        if (pthread_mutex_destroy(&(suffix##_##type##_MUTEX)) != 0) return 1;                           \
    }                                                                                                   \
    else if (suffix##_##type##_INIT == 0)                                                               \
    {                                                                                                   \
        suffix##_##type##_INIT = 1;                                                                     \
        suffix##_##type##_CALLBACK = fn;                                                                \
        if (pthread_mutex_init(&(suffix##_##type##_MUTEX), NULL) != 0) return 1;                        \
    }                                                                                                   \
    return 0;                                                                                           \
}                                                                                                       \
int _##suffix##_throw_##type(int level, int code, const char *file, int line, const char *fmt, ...)     \
{                                                                                                       \
    if(suffix##_##type##_INIT == 1 && suffix##_##type##_CALLBACK != NULL)                               \
    {                                                                                                   \
        if (pthread_mutex_lock(&(suffix##_##type##_MUTEX)) != 0) return 1;                              \
        va_list ap, ap_copy;                                                                            \
        va_start(ap, fmt);                                                                              \
        va_copy(ap_copy, ap);                                                                           \
        size_t len = vsnprintf(0, 0, fmt, ap_copy);                                                     \
        char *msg = malloc(len + 1);                                                                    \
        vsnprintf(msg, len + 1, fmt, ap);                                                               \
        (suffix##_##type##_CALLBACK)(level, code, file, line, msg);                                     \
        free(msg);                                                                                      \
        va_end(ap);                                                                                     \
        if (pthread_mutex_unlock(&(suffix##_##type##_MUTEX)) != 0) return 1;                            \
    }                                                                                                   \
    return 0;                                                                                           \
}                                                                                                       \

/*
*   LIST
*/

void _list_generic_create(void *list, size_t size, size_t elem_byte_size, const char *file, int line);
void _list_generic_destroy(void *list);
void _list_generic_resize(void *list, size_t elem_byte_size, size_t new_size);
void _list_generic_append(void *list, size_t elem_byte_size, const void *elem);
void _list_generic_remove(void *list, size_t elem_byte_size, size_t index);
void _list_generic_insert(void *list, size_t elem_byte_size, size_t index, const void *elem);
size_t _list_generic_find(void *list, size_t elem_byte_size, const void *elem);

#define LIST(type)                                                                                  \
struct List_##type                                                                                  \
{                                                                                                   \
    size_t size;                                                                                    \
    type *data;                                                                                     \
    size_t _capacity;                                                                               \
};                                                                                                  \
static inline struct List_##type _list_##type##_create(size_t size, const char *file, int line)     \
{                                                                                                   \
    struct List_##type list = {0};                                                                  \
    _list_generic_create(&list, size, sizeof(type), file, line);                                    \
    return list;                                                                                    \
}                                                                                                   \
static inline void list_##type##_destroy(struct List_##type *list)                                  \
{                                                                                                   \
    _list_generic_destroy(list);                                                                    \
}                                                                                                   \
static inline void list_##type##_resize(struct List_##type *list, size_t new_size)                  \
{                                                                                                   \
    _list_generic_resize(list, sizeof(type), new_size);                                             \
}                                                                                                   \
static inline void list_##type##_append(struct List_##type *list, const type *elem)                 \
{                                                                                                   \
    _list_generic_append(list, sizeof(type), elem);                                                 \
}                                                                                                   \
static inline void list_##type##_remove(struct List_##type *list, size_t index)                     \
{                                                                                                   \
    _list_generic_remove(list, sizeof(type), index);                                                \
}                                                                                                   \
static inline void list_##type##_insert(struct List_##type *list, size_t index, const type *elem)   \
{                                                                                                   \
    _list_generic_insert(list, sizeof(type), index, elem);                                          \
}                                                                                                   \
static inline size_t list_##type##_find(struct List_##type *list, const type *elem)                 \
{                                                                                                   \
    return _list_generic_find(list, sizeof(type), elem);                                            \
}                                                                                                   \
static inline size_t list_##type##_remove_on_finding(struct List_##type *list, const type *elem)    \
{                                                                                                   \
    size_t index = _list_generic_find(list, sizeof(type), elem);                                    \
    if (index != SIZE_MAX) _list_generic_remove(list, sizeof(type), index);                         \
    return index;                                                                                   \
}                                                                                                   \

#define list_create(type, size) _list_##type##_create(size, __FILE__, __LINE__)
#define for_list(index_name, list) for (size_t index_name = 0; index_name < list.size; ++index_name)

/*
* ARRAY
*/

#define Array(type, name, length, ...) struct { const size_t size; type data[length]; } name = { length, __VA_ARGS__ }

/*
*   MEM
*/

void *_mem_malloc(size_t size, const char *file, int line);
void *_mem_calloc(size_t nmemb, size_t size, const char *file, int line);
void *_mem_realloc(void *ptr, size_t size);
void _mem_free(void *ptr);

#define mem_malloc(size) _mem_malloc(size, __FILE__, __LINE__)
#define mem_calloc(nmemb, size) _mem_calloc(nmemb, size, __FILE__, __LINE__)
#define mem_realloc(ptr, size) _mem_realloc(ptr, size)
#define mem_free(ptr) _mem_free(ptr)
#define mem_alloc(type) mem_malloc(sizeof(type))
#define mem_alloc_init(type, name) type *name = mem_malloc(sizeof(type))

/*
*   TRACKER
*/

typedef void (*tracker_trace_fn)(void *address, const char *file, int line);

int tracker_register(void *address, const char *file, int line);
int tracker_unregister(void *address);
int tracker_change_register(void *old_address, void *new_address);
int tracker_trace(tracker_trace_fn fn);

#ifdef __cplusplus
}
#endif

#endif // __CUTILS_BASE_H__
