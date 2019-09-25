package de.unijena.bioinf.fingerid.connection_pooling;

import gnu.trove.procedure.TObjectProcedure;
import org.slf4j.LoggerFactory;

import java.io.Closeable;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.ReentrantLock;

/**
 * This class implements any kind of connection pool that manages some external resources (e.g. remote IO connections,
 * HTTP clients, SQL databases, ...).
 * <p>
 * The connection pool starts with a capacity that controls how many open connection are allowed to stay in the pool.
 * When first time a connection is requested, the connection pool will open a new connection and increase its size. When the connection is "closed", the connection pool will (instead of closing it) put the open connection into a list. When the next time a connection is requested, the connection pool will first empty its list before open a new connection. When the capacity of the pool is reached, the connection pool will block at every request until a free connection is available again.
 * <p>
 * Important note 1: The connection pool still requires that every requested connection is closed properly. The "real" closing of a connection is done when the connection pool itself is closed. For connections, that have a limited lifetime (e.g. SQL or HTTP connections), you might want to keep the connection pool in life only for a limited time, too, to avoid dying connection objects.
 * <p>
 * Important note 2: You should only order a single pooled connection per thread. Otherwise you might create a deadlock. The pool will print a warning in such cases
 * <blockquote><pre>
 * {@code
 * ConnectionPool<Resource> pool = new ConnectionPool(resourceConnector, 10);
 * for (int i=0; i < 1000; ++i) {
 *   // even if 1000 threads order connections, there will be never more than 10 connection objects
 *   // ever created
 *   new Thread(new Runnable() {
 *       try (final Resource r = pool.orderConnection()) {
 *         // dome something with resource
 *         // inside this block, the resource is exclusively available for this thread. No other thread
 *         // can access this resource
 *       }
 *   }).start();
 * }
 * pool.close(); // will close all resources
 * }
 * </pre></blockquote>
 *
 * @param <T>
 */
public class ConnectionPool<T> implements Closeable, AutoCloseable {

    public interface Connection<T> {
        public T open() throws IOException;

        public void close(T connection) throws IOException;
    }

    protected final Connection<T> connector;
    protected final ConcurrentLinkedQueue<T> freeConnections;
    protected volatile boolean shutdown, forcedShutdown;
    protected final AtomicInteger size, sharedCounter;
    protected final int capacity;
    protected final Condition noFreeConnectionsLeft, noOpenConnections;
    protected final ReentrantLock connectionLock;
    protected final AtomicInteger waitingThreads;

    protected volatile long lastConnectionCheck;

    private final Set<Thread> threads = Collections.newSetFromMap(new ConcurrentHashMap<>());

    public ConnectionPool(Connection<T> connector, int capacity) {
        this.connector = connector;
        this.freeConnections = new ConcurrentLinkedQueue<>();
        this.shutdown = false;
        this.size = new AtomicInteger(0);
        this.capacity = capacity;
        this.connectionLock = new ReentrantLock();
        this.noFreeConnectionsLeft = connectionLock.newCondition();
        this.noOpenConnections = connectionLock.newCondition();
        this.waitingThreads = new AtomicInteger(0);
        this.sharedCounter = new AtomicInteger(1);
        this.lastConnectionCheck = System.currentTimeMillis();
    }

    /**
     * if since the last time this method was called numberOfSeconds seconds are passed, close all idling connections.
     */
    public void testConnectionsAfter(int numberOfSeconds) throws IOException {
        final long time = System.currentTimeMillis();
        final long check = time - numberOfSeconds;
        if (lastConnectionCheck < check) {
            synchronized (this) { // prevent that closeAllIdlingConnections is called by many threads in parallel
                if (lastConnectionCheck < time) {
                    closeAllIdlingConnections();
                    lastConnectionCheck = System.currentTimeMillis();
                }
            }
        }
    }

    protected int getSize() {
        return size.get();
    }

    public ConnectionPool(Connection<T> connector) {
        this(connector, Integer.MAX_VALUE);
    }

    public PooledConnection<T> orderConnection() throws InterruptedException, IOException {

        if (shutdown)
            throw new IllegalStateException("Connection pool is closed and does not accept new requests.");
        if (threads.contains(Thread.currentThread()))
            LoggerFactory.getLogger(getClass()).warn(Thread.currentThread().getName() + " has already a connection of. " + toString() + " Ordering multiple connections with the same thread might cause a deadlock. See Stacktrace:" + System.lineSeparator() + Arrays.toString(Thread.currentThread().getStackTrace()));
        waitingThreads.incrementAndGet();
        final PooledConnection<T> c = orderConnectionDontIncreaseWaitingCount();
        threads.add(Thread.currentThread());
        return c;

    }

    private PooledConnection<T> orderConnectionDontIncreaseWaitingCount() throws InterruptedException, IOException {
        final T connection = freeConnections.poll();
        if (connection == null) {
            // try to open a new connection
            final boolean nuConnection;
            synchronized (size) {
                nuConnection = size.get() < capacity;
                if (nuConnection)
                    size.incrementAndGet();
            }

            if (nuConnection) {//return new connection added to the pool
                try {
                    return new PooledConnection<T>(this, connector.open());
                } catch (IOException e) {
                    connectionLock.lock();
                    size.decrementAndGet();
                    try {
                        noFreeConnectionsLeft.signalAll(); // we might get new capacity free
                    } finally {
                        connectionLock.unlock();
                    }
                    throw e;
                }
            } else {//wait for connection
                return waitForNewConnectionComesIn();
            }
        } else {//return free existing connection
            return new PooledConnection<>(this, connection);
        }
    }


    /**
     * This will close all idling connections. This method might be called regularly for e.g. database connections
     * to avoid that a database connection (that is idling for a long time) is dying.
     */
    public void closeAllIdlingConnections() throws IOException {
        while (!freeConnections.isEmpty()) {
            final T c = freeConnections.poll();
            if (c != null) {
                connector.close(c);
                size.decrementAndGet();
            }
        }
        connectionLock.lock();
        try {
            noFreeConnectionsLeft.signalAll(); // we might get new capacity free
        } finally {
            connectionLock.unlock();
        }
    }

    /**
     * This will refresh all idling connections. This method might be called regularly for e.g. database connections
     * to avoid that a database connection (that is idling for a long time) is dying.
     */
    public void refreshAllIdlingConnections(TObjectProcedure<T> refreshOperation) {
        final ArrayList<T> refreshedConnections = new ArrayList<>();
        while (!freeConnections.isEmpty()) {
            final T c = freeConnections.poll();
            if (c != null) {
                refreshOperation.execute(c);
                refreshedConnections.add(c);
            }
        }
        freeConnections.addAll(refreshedConnections);
        connectionLock.lock();
        try {
            noFreeConnectionsLeft.signalAll(); // we might get new capacity free
        } finally {
            connectionLock.unlock();
        }
    }

    private PooledConnection<T> waitForNewConnectionComesIn() throws InterruptedException, IOException {
//        System.out.println(Thread.currentThread().getName() + " waits for connection of " + connector.getClass().getName());
        while (true) {
            connectionLock.lock();
            try {
                noFreeConnectionsLeft.await();
                if (forcedShutdown) {
//                    System.out.println(Thread.currentThread().getName() + ": waiting for connection DONE (with Exception)!");
                    throw new InterruptedException("Interrupted by shutdown of connection pool");
                }
                final T c = freeConnections.poll();
                if (c != null) {
//                    System.out.println(Thread.currentThread().getName() + ": waiting for connection DONE (free Conn. available)!");
                    return new PooledConnection<T>(this, c); //this refreshs a null connection
                }
            } finally {
                connectionLock.unlock();
            }

            if (size.get() < capacity) {
//                System.out.println(Thread.currentThread().getName() + ": waiting for connection DONE (SIZE < CAP)!");
                return orderConnectionDontIncreaseWaitingCount();
            }
        }
    }

    void freeConnection(final PooledConnection<T> connection) throws IOException {
        try {
//            System.out.println("---> " + Thread.currentThread().getName() + " START freeing " + connection.connection.toString() + " (size: " + size.get() + "/" + capacity + ")");
            if (connection == null) throw new NullPointerException();
            if (connection.closed) return; // already freed
            synchronized (connection) {
                if (connection.closed) return;
                connection.closed = true;
            }

            if (waitingThreads.decrementAndGet() <= 0) {
                connectionLock.lock();
                noOpenConnections.signalAll();
                connectionLock.unlock();
            }
            // else move connection into free list
            if (forcedShutdown) {
                connector.close(connection.connection);
                return;
            }
            freeConnections.add(connection.connection);
            connectionLock.lock();
            noFreeConnectionsLeft.signal();
            connectionLock.unlock();
//            System.out.println("---> " + Thread.currentThread().getName() + " END freeing " + connection.connection.toString() + " (size: " + size.get() + "/" + capacity + ")");
        } finally {
            threads.remove(Thread.currentThread());
        }
    }

    public void shutdown() throws InterruptedException, IOException {
        if (sharedCounter.decrementAndGet() > 0) return;
        shutdown = true;
        while (true) {
            connectionLock.lock();
            if (waitingThreads.get() <= 0) {
                connectionLock.unlock();
                break;
            }
            noOpenConnections.await();
            connectionLock.unlock();
        }
        // close all connections
        while (!freeConnections.isEmpty()) {
            connector.close(freeConnections.poll());
        }
    }

    @Override
    public void close() throws IOException {
        // this method blocks until all requests are processed
        try {
            shutdown();
        } catch (InterruptedException e) {
            enforceShutdown();
        }
    }

    /**
     * This will terminate all running connections. However, it cannot
     * guarantee that all open connections are properly closed.
     */
    public void enforceShutdown() throws IOException {
        if (sharedCounter.decrementAndGet() > 0) return;
        shutdown = true;
        forcedShutdown = true;
        while (!freeConnections.isEmpty()) {
            connector.close(freeConnections.poll());
        }
        connectionLock.lock();
        try {
            noFreeConnectionsLeft.signalAll(); // we might get new capacity free
        } finally {
            connectionLock.unlock();
        }
    }

    public ConnectionPool<T> newSharedConnectionPool() {
        sharedCounter.incrementAndGet();
        return this;
    }

    public int getCapacity() {
        return capacity;
    }

    public int getNumberOfIdlingConnections() {
        return freeConnections.size();
    }

    @Override
    public String toString() {
        return super.toString() + "{" + connector.getClass().getName() + "}";
    }
}
