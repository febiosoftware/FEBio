// #include <stdio.h> 
// #include <string.h>
// #include <thread>
// #ifdef WIN32
// #include <winsock2.h>
// #include <ws2tcpip.h>
// #else
// #include <sys/types.h> 
// #include <unistd.h> 
// #include <sys/socket.h> 
// #include <arpa/inet.h> 
// #include <netinet/in.h> 
// #include <netdb.h>
// #endif
#include "ping.h"
// #include <FEBioLib/version.h>

// #define PORT     8560 
// #define MSGLEN 50 

// #ifdef WIN32
// #define PLAT "Windows"
// #elif __APPLE__
// #define PLAT "macOS"
// #else
// #define PLAT "Linux"
// #endif

// // Main work done in a thread so that the DNS lookup doesn't pause execution 
// void pingThread() {
// 	int sockfd;
// 	struct sockaddr_in servaddr;

//     char message[MSGLEN];
//     sprintf(message, "%s %d.%d.%d", PLAT, VERSION, SUBVERSION, SUBSUBVERSION);

// #ifdef WIN32
// 	WSADATA wsaData;

// 	int iResult;

// 	// Initialize Winsock
// 	iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
// 	if (iResult != 0) {
// 		return;
// 	}
// #endif

// 	// Creating socket descriptor 
// 	sockfd = socket(AF_INET, SOCK_DGRAM, 0);
// #ifdef WIN32
// 	if(sockfd == INVALID_SOCKET)
// #else
// 	if (sockfd  < 0) 
// #endif
// 	{
// 		return;
// 	}

// 	// Strucs for DNS lookup
// 	struct addrinfo *result = NULL;
// 	struct addrinfo hints;

// 	//--------------------------------
// 	// Setup the hints address info structure
// 	// which is passed to the getaddrinfo() function
// 	memset(&hints, 0, sizeof(hints));
// 	hints.ai_family = AF_UNSPEC;
// 	hints.ai_socktype = SOCK_STREAM;
// 	hints.ai_protocol = IPPROTO_TCP;

// 	const char address[] = "repo.febio.org";
// 	const char port[] = "";

// 	// DNS Lookup
// 	int dwRetval = getaddrinfo(address, port, &hints, &result);
// 	if (dwRetval == 0 && result->ai_family == AF_INET) 
//     {
//         // Filling server information 
//         memset(&servaddr, 0, sizeof(servaddr));
//         servaddr.sin_family = AF_INET;
//         servaddr.sin_port = htons(PORT);
//         servaddr.sin_addr.s_addr = ((struct sockaddr_in *) result->ai_addr)->sin_addr.s_addr;

//         sendto(sockfd, (const char *)message, strlen(message),
//             0, (const struct sockaddr *) &servaddr,
//             sizeof(servaddr));

//         freeaddrinfo(result);
// 	}

// #ifdef WIN32
// 	closesocket(sockfd);
//     WSACleanup();
// #else
// 	close(sockfd);
// #endif
// }

void ping()
{
    // std::thread (pingThread).detach();
}